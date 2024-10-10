/***********************************************************************
 * Copyright (c) 2024 Jayamine Alupotha                                  *
 * Distributed under the MIT software license, see the accompanying    *
 * file COPYING2 or http://www.opensource.org/licenses/mit-license.php. *
 ***********************************************************************/

#ifndef RAHAS_SECP256K1_MAIN_IMPL_H
#define RAHAS_SECP256K1_MAIN_IMPL_H

#include <openssl/rand.h>
#include "group.h"
#include "modules/commitment/main_impl.h"
#include "math.h"
#include "testrand_impl.h"
#include "modules/ringcip/main_impl.h"
#include "main_impl.h"
#include "zero_com_impl.h"

static void secp256k1_ge_save(unsigned char *buf, secp256k1_ge* ge) {
    secp256k1_fe_normalize(&ge->x);
    secp256k1_fe_get_b32(&buf[1], &ge->x);
    buf[0] = 9 ^ secp256k1_fe_is_quad_var(&ge->y);
}

static void secp256k1_ge_load(unsigned char *buf, secp256k1_ge* ge) {
    secp256k1_fe fe;
    secp256k1_fe_set_b32(&fe, &buf[1]);
    secp256k1_ge_set_xquad(ge, &fe);
    if (buf[0] & 1) {
        secp256k1_ge_neg(ge, ge);
    }
}

int secp256k1_identity_create(
        const secp256k1_context* ctx,
        secp256k1_pedersen_commitment *com,
        const unsigned char *blind,
        secp256k1_hash digest,
        const secp256k1_generator *value_gen,
        const secp256k1_generator *blind_gen
) {
    return secp256k1_pedersen_blind_commit(ctx, com, blind, digest.data, value_gen, blind_gen);
}


int secp256k1_get_multi_generators(
        const secp256k1_context* ctx,
        secp256k1_generator *gens,
        const unsigned char *gen_seed,
        int m
) {
    secp256k1_sha256 sha;
    unsigned char buf[32];
    int i;

    memset(buf, 0, 32);
    for (i = 0; i < m; i++) {
        secp256k1_sha256_initialize(&sha);
        secp256k1_sha256_write(&sha, gen_seed, 32);
        secp256k1_sha256_write(&sha, buf, 32);
        secp256k1_sha256_finalize(&sha, buf);

        if(secp256k1_generator_generate(ctx, &gens[i], buf) == 0)
            i--;
    }
    return 1;
}


int secp256k1_get_multi_gen_commitment(
        secp256k1_pedersen_commitment *com,
        const unsigned char *blind,
        const unsigned char *digests,
        const secp256k1_generator *blind_gen,
        const secp256k1_generator *gens,
        int m
) {
    secp256k1_scalar tmp;
    secp256k1_ge tmpG;
    secp256k1_gej tmpGj;
    secp256k1_gej tmpGj1;
    int i, overflow;

    secp256k1_generator_load(&tmpG, blind_gen);
    secp256k1_scalar_set_b32(&tmp, blind, &overflow);
    if (overflow != 0)
        return 0;
    secp256k1_ecmult_const(&tmpGj, &tmpG, &tmp, 256);
    for (i = 0; i < m; i++) {
        secp256k1_generator_load(&tmpG, &gens[i]);
        secp256k1_scalar_set_b32(&tmp, digests + i * 32, &overflow);
        if (secp256k1_scalar_is_zero(&tmp))
            continue;
        secp256k1_ecmult_const(&tmpGj1, &tmpG, &tmp, 256);
        secp256k1_ge_set_gej(&tmpG, &tmpGj1);
        secp256k1_gej_add_ge(&tmpGj, &tmpGj, &tmpG);
    }
    secp256k1_ge_set_gej(&tmpG, &tmpGj);
    secp256k1_fe_normalize(&tmpG.x);
    secp256k1_pedersen_commitment_save(com, &tmpG);

    return 1;
}


int secp256k1_unlinked_logarithmic_zero_com_prove(
        const secp256k1_context* ctx,
        unsigned char *proof,
        int index,
        const unsigned char *blind,
        int ring_size,
        int n,
        const unsigned char *challenge,
        const secp256k1_pedersen_commitment *coms,
        const secp256k1_generator *value_gen,
        const secp256k1_generator *blind_gen
) {
    secp256k1_pedersen_commitment J[n] ;
    secp256k1_pedersen_commitment A[n];
    secp256k1_pedersen_commitment B[n];
    secp256k1_pedersen_commitment D[n];
    secp256k1_scalar f[n];
    secp256k1_scalar za[n];
    secp256k1_scalar zb[n];
    secp256k1_scalar zd;

    int N = (1 << n);

    secp256k1_scalar r[n];
    secp256k1_scalar a[n];
    secp256k1_scalar s[n];
    secp256k1_scalar v[n];
    secp256k1_scalar rho[n];
    uint8_t r32[n][32];
    uint8_t a32[n][32];
    uint8_t s32[n][32];
    uint8_t v32[n][32];
    uint8_t rho32[n][32];
    secp256k1_scalar d[n][2];
    secp256k1_scalar p[N][n + 1];
    secp256k1_scalar tmp1;
    secp256k1_scalar tmp2;
    secp256k1_scalar tmp3;
    secp256k1_ge geng;
    secp256k1_ge genh;
    secp256k1_gej Cj;
    secp256k1_gej Dj;
    secp256k1_ge tmpG;
    secp256k1_scalar x32;
    secp256k1_sha256 sha;
    unsigned char buf[32];
    int i, l, l1, overflow, i_l, j_l;
    int pointer = 0;

    VERIFY_CHECK(ctx != NULL);
    ARG_CHECK(secp256k1_ecmult_context_is_built(&ctx->ecmult_ctx));
    ARG_CHECK(proof != NULL);
    ARG_CHECK(challenge != NULL);
    ARG_CHECK(ring_size != 0);
    ARG_CHECK(n != 0);
    ARG_CHECK(!(index > N || index < 0));

    secp256k1_generator_load(&geng, blind_gen);
    secp256k1_generator_load(&genh, value_gen);

    /* Set random challenges */
    for (l = 0; l < n; l++) {
        secp256k1_rand256(r32[l]);
        secp256k1_scalar_set_b32(&r[l], r32[l], &overflow);
        secp256k1_rand256(a32[l]);
        secp256k1_scalar_set_b32(&a[l], a32[l], &overflow);
        secp256k1_rand256(s32[l]);
        secp256k1_scalar_set_b32(&s[l], s32[l], &overflow);
        secp256k1_rand256(v32[l]);
        secp256k1_scalar_set_b32(&v[l], v32[l], &overflow);
        secp256k1_rand256(rho32[l]);
        secp256k1_scalar_set_b32(&rho[l], rho32[l], &overflow);
    }

    /* Find polynomial coefficients */
    /*
     * P(i) = prod F(index, i_j) = prod d(i_j, l_j)Z + (-1)^{d_{0, i_j}}a(index)
     */
    for (i = 0; i < N; i++) {
        /* Set p[i][l] = 0 */
        for (l = 0; l < n + 1; l++) {
            secp256k1_scalar_set_int(&p[i][l], 0);
        }
        for (l = 0; l < n; l++) {
            secp256k1_scalar_set_int(&d[l][0], 0);
            secp256k1_scalar_set_int(&d[l][1], 0);
        }

        /* Multiply by d(i_j, l_j)Z + (-1)^{d_{0, i_j}}a(index) */
        for (l = 0; l < n; l++) {
            i_l = (i & (1 << l)) >> l;
            j_l = (index & (1 << l)) >> l;

            secp256k1_scalar_set_b32(&d[l][0], a32[l], &overflow);
            if (i_l == 0) {
                secp256k1_scalar_negate(&d[l][0], &d[l][0]);
                secp256k1_scalar_set_int(&d[l][1], 1 - j_l);
            } else {
                secp256k1_scalar_set_int(&d[l][1], j_l);
            }
        }

        secp256k1_scalar_add(&p[i][0], &p[i][0], &d[0][0]);
        secp256k1_scalar_add(&p[i][1], &p[i][1], &d[0][1]);
        for (l = 1; l < n; l++) {
            secp256k1_scalar_set_int(&tmp1, 0);
            secp256k1_scalar_add(&tmp1, &tmp1, &p[i][0]);
            secp256k1_scalar_mul(&p[i][0], &p[i][0], &d[l][0]);
            for (l1 = 0; l1 < l; l1++) {
                /* tmp3 = d[l][1] * p[i][l] (X^(l + 1)) */
                secp256k1_scalar_mul(&tmp1, &tmp1, &d[l][1]);
                /* tmp4 = d[l][0] * p[i][l + 1] (X^(l + 1)) */
                secp256k1_scalar_mul(&tmp2, &p[i][l1 + 1], &d[l][0]);
                secp256k1_scalar_add(&tmp2, &tmp1, &tmp2);
                /* Copy &p[i][l1 + 1] */
                secp256k1_scalar_set_int(&tmp1, 0);
                secp256k1_scalar_add(&tmp1, &tmp1, &p[i][l1 + 1]);
                /* Update &p[i][l + 1]*/
                secp256k1_scalar_set_int(&p[i][l1 + 1], 0);
                secp256k1_scalar_add(&p[i][l1 + 1], &p[i][l1 + 1], &tmp2);
            }
            secp256k1_scalar_mul(&p[i][l + 1], &tmp1, &d[l][1]);
        }
    }

    /*
     * J_l  = com(j_l, r_l)
     * A_l  = com(a_l, s_l)
     * B_l  = com(a_l * j_l, t_l)
     * D_l = prod_{i=1}^{N} (𝐶𝑖ℎ−𝑑 )𝑝𝑖−1,𝑙 −1)× COM.cmt(0, 𝜌𝑙−1)
     */
    for (l = 0; l < n; l++) {
        j_l = (index & (1 << l)) >> l;
        if (!secp256k1_pedersen_commit(ctx, &J[l], r32[l], j_l,
                                       value_gen,
                                       blind_gen)
            || !secp256k1_pedersen_blind_commit(ctx, &A[l], s32[l], a32[l],
                                                value_gen,
                                                blind_gen)) {
            return 0;
        }

        if (j_l == 1) {
            if (!secp256k1_pedersen_blind_commit(ctx, &B[l], v32[l], a32[l],
                                                 value_gen,
                                                 blind_gen)) {
                return 0;
            }
        } else {
            if (!secp256k1_pedersen_commit(ctx, &B[l], v32[l], 0,
                                           value_gen,
                                           blind_gen)) {
                return 0;
            }
        }

        if (!secp256k1_pedersen_commit(ctx, &D[l], rho32[l], 0,
                                       value_gen,
                                       blind_gen)) {
            return 0;
        }

        secp256k1_pedersen_commitment_load(&tmpG, &D[l]);
        secp256k1_gej_set_ge(&Dj, &tmpG);
        for (i = 0; i < ring_size; i++) {
            if (secp256k1_scalar_is_zero(&p[i][l]) == 1) {
                continue;
            }
            secp256k1_pedersen_commitment_load(&tmpG, &coms[i]);
            secp256k1_ecmult_const(&Cj, &tmpG, &p[i][l], 256);
            secp256k1_ge_set_gej(&tmpG, &Cj);
            secp256k1_gej_add_ge(&Dj, &Dj, &tmpG);
        }
        secp256k1_ge_set_gej(&tmpG, &Dj);
        secp256k1_fe_normalize(&tmpG.x);
        secp256k1_pedersen_commitment_save(&D[l], &tmpG);

        secp256k1_pedersen_commitment_serialize(ctx, proof + pointer, &J[l]);
        pointer += 33;
        secp256k1_pedersen_commitment_serialize(ctx, proof + pointer, &A[l]);
        pointer += 33;
        secp256k1_pedersen_commitment_serialize(ctx, proof + pointer, &B[l]);
        pointer += 33;
        secp256k1_pedersen_commitment_serialize(ctx, proof + pointer, &D[l]);
        pointer += 33;

    }

    /* x */
    secp256k1_sha256_initialize(&sha);
    secp256k1_sha256_write(&sha, challenge, 32);
    secp256k1_sha256_write(&sha, proof, 33 * n * 4);
    secp256k1_sha256_finalize(&sha, buf);
    secp256k1_scalar_set_b32(&x32, buf, NULL);


    secp256k1_scalar_set_int(&tmp1, 1);
    secp256k1_scalar_set_int(&tmp2, 0);
    for (l = 0; l < n; l++) {
        j_l = (index & (1 << l)) >> l;
        secp256k1_scalar_set_b32(&f[l], a32[l], &overflow);
        if (j_l == 1) {
            secp256k1_scalar_add(&f[l], &f[l], &x32);
        }

        secp256k1_scalar_set_b32(&za[l], buf, &overflow);
        secp256k1_scalar_negate(&zb[l], &f[l]);
        secp256k1_scalar_add(&zb[l], &zb[l], &x32);

        secp256k1_scalar_mul(&za[l], &za[l], &r[l]);
        secp256k1_scalar_mul(&zb[l], &zb[l], &r[l]);

        secp256k1_scalar_add(&za[l], &za[l], &s[l]);
        secp256k1_scalar_add(&zb[l], &zb[l], &v[l]);

        secp256k1_scalar_mul(&tmp3, &rho[l], &tmp1);
        secp256k1_scalar_add(&tmp2, &tmp2, &tmp3);

        secp256k1_scalar_mul(&tmp1, &tmp1, &x32); /* should be after tmp2 */

        secp256k1_scalar_get_b32(proof + pointer, &f[l]);
        pointer += 32;
        secp256k1_scalar_get_b32(proof + pointer, &za[l]);
        pointer += 32;
        secp256k1_scalar_get_b32(proof + pointer, &zb[l]);
        pointer += 32;
    }
    secp256k1_scalar_set_b32(&zd, blind, &overflow); /* the known secret key */
    secp256k1_scalar_mul(&zd, &zd, &tmp1);
    secp256k1_scalar_negate(&tmp2, &tmp2);
    secp256k1_scalar_add(&zd, &zd, &tmp2);
    secp256k1_scalar_get_b32(proof + pointer, &zd);


    for (l = 0; l < n; l++) {
        secp256k1_rand256(a32[l]);
        memset(r32[l], 0, 32);
        memset(r32[l], 0, 32);
        memset(s32[l], 0, 32);
        memset(v32[l], 0, 32);
        memset(rho32[l], 0, 32);
        secp256k1_scalar_set_int(&a[l], 0);
        secp256k1_scalar_set_int(&r[l], 0);
        secp256k1_scalar_set_int(&s[l], 0);
        secp256k1_scalar_set_int(&v[l], 0);
        secp256k1_scalar_set_int(&rho[l], 0);
    }

    return 1;
}

int secp256k1_unlinked_logarithmic_zero_com_verify(
        const secp256k1_context* ctx,
        unsigned char *proof,
        int ring_size,
        int n,
        const unsigned char *challenge,
        const secp256k1_pedersen_commitment *coms,
        const secp256k1_generator *value_gen,
        const secp256k1_generator *blind_gen
) {
    secp256k1_pedersen_commitment J[n];
    secp256k1_pedersen_commitment A[n];
    secp256k1_pedersen_commitment B[n];
    secp256k1_pedersen_commitment D[n];
    secp256k1_scalar f[n];
    secp256k1_scalar za[n];
    secp256k1_scalar zb[n];
    secp256k1_scalar zd;

    int N = (1 << n);

    secp256k1_scalar tmpf[N];
    secp256k1_scalar tmp2;
    secp256k1_scalar tmp3;
    secp256k1_scalar tmp4;
    secp256k1_ge geng;
    secp256k1_ge genh;
    secp256k1_gej tmpj1;
    secp256k1_gej tmpj2;
    secp256k1_gej tmpj3;
    secp256k1_ge tmpG1;
    secp256k1_ge tmpG2;
    secp256k1_pedersen_commitment com;
    secp256k1_scalar x32;
    secp256k1_sha256 sha;
    unsigned char buf[32];
    unsigned char buf_blind[32];
    unsigned char RHS[32];
    unsigned char LHS[32];
    int i, l, overflow, i_l;
    int pointer = 0;

    VERIFY_CHECK(ctx != NULL);
    ARG_CHECK(secp256k1_ecmult_context_is_built(&ctx->ecmult_ctx));
    ARG_CHECK(proof != NULL);
    ARG_CHECK(challenge != NULL);
    ARG_CHECK(N != 0);
    ARG_CHECK(n != 0);

    secp256k1_generator_load(&geng, blind_gen);
    secp256k1_generator_load(&genh, value_gen);

    /* x */
    secp256k1_sha256_initialize(&sha);
    secp256k1_sha256_write(&sha, challenge, 32);
    secp256k1_sha256_write(&sha, proof, 33 * n * 4);
    secp256k1_sha256_finalize(&sha, buf);
    secp256k1_scalar_set_b32(&x32, buf, NULL);

    for (l = 0; l < n; l++) {
        if (secp256k1_pedersen_commitment_parse(ctx, &J[l], proof + pointer) == 0)
            return 0;
        pointer += 33;
        if (secp256k1_pedersen_commitment_parse(ctx, &A[l], proof + pointer) == 0)
            return 0;
        pointer += 33;
        if (secp256k1_pedersen_commitment_parse(ctx, &B[l], proof + pointer) == 0)
            return 0;
        pointer += 33;
        if (secp256k1_pedersen_commitment_parse(ctx, &D[l], proof + pointer) == 0)
            return 0;
        pointer += 33;
    }

    for (l = 0; l < n; l++) {
        secp256k1_scalar_set_b32(&f[l], proof + pointer, &overflow);
        pointer += 32;
        secp256k1_scalar_set_b32(&za[l], proof + pointer, &overflow);
        pointer += 32;
        secp256k1_scalar_set_b32(&zb[l], proof + pointer, &overflow);
        pointer += 32;
    }
    secp256k1_scalar_set_b32(&zd, proof + pointer, &overflow);

    /* Verification */
    for (l = 0; l < n; l++) {
        secp256k1_pedersen_commitment_load(&tmpG1, &J[l]);
        secp256k1_ecmult_const(&tmpj1, &tmpG1, &x32, 256);
        secp256k1_pedersen_commitment_load(&tmpG2, &A[l]);
        secp256k1_gej_add_ge(&tmpj1, &tmpj1, &tmpG2);

        secp256k1_scalar_get_b32(buf, &f[l]);
        secp256k1_scalar_get_b32(buf_blind, &za[l]);
        if (!secp256k1_pedersen_blind_commit(ctx, &com, buf_blind, buf,
                                             value_gen,
                                             blind_gen)) {
            return 0;
        }
        secp256k1_ge_set_gej(&tmpG1, &tmpj1);
        secp256k1_fe_normalize(&tmpG1.x);
        secp256k1_fe_get_b32(LHS, &tmpG1.x);

        secp256k1_pedersen_commitment_load(&tmpG2, &com);
        secp256k1_fe_normalize(&tmpG2.x);
        secp256k1_fe_get_b32(RHS, &tmpG2.x);

        if (memcmp(LHS, RHS, 32) != 0)
            return 0;

        secp256k1_pedersen_commitment_load(&tmpG1, &J[l]);
        secp256k1_scalar_negate(&tmp3, &f[l]);
        secp256k1_scalar_add(&tmp3, &tmp3, &x32);
        secp256k1_ecmult_const(&tmpj1, &tmpG1, &tmp3, 256);
        secp256k1_pedersen_commitment_load(&tmpG2, &B[l]);
        secp256k1_gej_add_ge(&tmpj1, &tmpj1, &tmpG2);

        secp256k1_scalar_get_b32(buf_blind, &zb[l]);
        if (!secp256k1_pedersen_commit(ctx, &com, buf_blind, 0,
                                       value_gen,
                                       blind_gen)) {
            return 0;
        }
        secp256k1_ge_set_gej(&tmpG1, &tmpj1);
        secp256k1_fe_normalize(&tmpG1.x);
        secp256k1_fe_get_b32(LHS, &tmpG1.x);

        secp256k1_pedersen_commitment_load(&tmpG2, &com);
        secp256k1_fe_normalize(&tmpG2.x);
        secp256k1_fe_get_b32(RHS, &tmpG2.x);

        if (memcmp(LHS, RHS, 32) != 0)
            return 0;

    }


    /*
     * Correctness of the polynomial coefficients
     */
    for (i = 0; i < N; i++) {
        secp256k1_scalar_set_int(&tmpf[i], 1);
        for (l = 0; l < n; l++) {
            i_l = (i & (1 << l)) >> l;

            if (i_l == 1) {
                secp256k1_scalar_mul(&tmpf[i], &tmpf[i], &f[l]);
            } else {
                secp256k1_scalar_negate(&tmp4, &f[l]);
                secp256k1_scalar_add(&tmp4, &tmp4, &x32);
                secp256k1_scalar_mul(&tmpf[i], &tmpf[i], &tmp4);
            }
        }
    }

    /* ----------------------------------------------------------------------- */

    secp256k1_gej_set_infinity(&tmpj3);
    for (i = 0; i < ring_size; i++) {

        secp256k1_pedersen_commitment_load(&tmpG1, &coms[i]);
        secp256k1_ecmult_const(&tmpj2, &tmpG1, &tmpf[i], 256);
        secp256k1_ge_set_gej(&tmpG1, &tmpj2);
        secp256k1_gej_add_ge(&tmpj3, &tmpj3, &tmpG1);
    }

    secp256k1_scalar_set_int(&tmp2, 1);
    for (l = 0; l < n; l++) {
        secp256k1_scalar_negate(&tmp3, &tmp2);
        secp256k1_pedersen_commitment_load(&tmpG2, &D[l]);
        secp256k1_ecmult_const(&tmpj2, &tmpG2, &tmp3, 256);
        secp256k1_ge_set_gej(&tmpG2, &tmpj2);
        secp256k1_gej_add_ge(&tmpj3, &tmpj3, &tmpG2);

        secp256k1_scalar_mul(&tmp2, &tmp2, &x32);
    }

    secp256k1_scalar_get_b32(buf_blind, &zd);
    if (!secp256k1_pedersen_commit(ctx, &com, buf_blind, 0,
                                   value_gen,
                                   blind_gen)) {
        return 0;
    }
    secp256k1_ge_set_gej(&tmpG1, &tmpj3);
    secp256k1_fe_normalize(&tmpG1.x);
    secp256k1_fe_get_b32(RHS, &tmpG1.x);

    secp256k1_pedersen_commitment_load(&tmpG2, &com);
    secp256k1_fe_normalize(&tmpG2.x);
    secp256k1_fe_get_b32(LHS, &tmpG2.x);

    if (memcmp(LHS, RHS, 32) != 0)
        return 0;

    return 1;
}

int secp256k1_unlinked_logarithmic_identity_prove(
        const secp256k1_context* ctx,
        secp256k1_unlinked_identity_pf *proof,
        int index,
        const unsigned char *blind,
        const unsigned char *digest,
        const unsigned char *challenge,
        int ring_size,
        int n,
        secp256k1_pedersen_commitment *coms,
        const secp256k1_generator *value_gen,
        const secp256k1_generator *blind_gen
) {
    secp256k1_pedersen_commitment C[ring_size];

    int N = (1 << n);

    secp256k1_ge geng;
    secp256k1_ge genh;
    secp256k1_gej hnegj;
    secp256k1_gej Cj;
    secp256k1_ge tmpG;
    secp256k1_scalar negdigest32;
    int i, overflow;

    VERIFY_CHECK(ctx != NULL);
    ARG_CHECK(secp256k1_ecmult_context_is_built(&ctx->ecmult_ctx));
    ARG_CHECK(proof != NULL);
    ARG_CHECK(digest != NULL);
    ARG_CHECK(challenge != NULL);
    ARG_CHECK(ring_size != 0);
    ARG_CHECK(n != 0);
    ARG_CHECK(!(index > N || index < 0));

    secp256k1_generator_load(&geng, blind_gen);
    secp256k1_generator_load(&genh, value_gen);

    secp256k1_scalar_set_b32(&negdigest32, digest, &overflow);
    secp256k1_scalar_negate(&negdigest32, &negdigest32);
    secp256k1_ecmult_const(&hnegj, &genh, &negdigest32, 256);

    proof->data = (uint8_t *) malloc((32 * (n * 3 + 1) + 33 * n * 4) * sizeof(uint8_t));

    for (i = 0; i < ring_size; i++) {
        secp256k1_pedersen_commitment_load(&tmpG, &coms[i]);
        secp256k1_gej_add_ge(&Cj, &hnegj, &tmpG);
        secp256k1_ge_set_gej(&tmpG, &Cj);
        secp256k1_pedersen_commitment_save(&C[i], &tmpG);
    }

    return secp256k1_unlinked_logarithmic_zero_com_prove(ctx, proof->data, index, blind,
                                                         ring_size, n, challenge, C, value_gen, blind_gen);
}

int secp256k1_unlinked_logarithmic_identity_verify(
        const secp256k1_context* ctx,
        secp256k1_unlinked_identity_pf *proof,
        const unsigned char *digest,
        const unsigned char *challenge,
        secp256k1_pedersen_commitment *coms,
        int ring_size,
        int n,
        const secp256k1_generator *value_gen,
        const secp256k1_generator *blind_gen
) {
    secp256k1_pedersen_commitment C[ring_size];

    int N = (1 << n);

    secp256k1_ge geng;
    secp256k1_ge genh;
    secp256k1_gej hnegj;
    secp256k1_gej Cj;
    secp256k1_ge tmpG;
    secp256k1_scalar negdigest32;
    int i, overflow;

    VERIFY_CHECK(ctx != NULL);
    ARG_CHECK(secp256k1_ecmult_context_is_built(&ctx->ecmult_ctx));
    ARG_CHECK(proof != NULL);
    ARG_CHECK(digest != NULL);
    ARG_CHECK(challenge != NULL);
    ARG_CHECK(N != 0);
    ARG_CHECK(n != 0);

    secp256k1_generator_load(&geng, blind_gen);
    secp256k1_generator_load(&genh, value_gen);

    secp256k1_scalar_set_b32(&negdigest32, digest, &overflow);
    secp256k1_scalar_negate(&negdigest32, &negdigest32);
    secp256k1_ecmult_const(&hnegj, &genh, &negdigest32, 256);

    for (i = 0; i < ring_size; i++) {
        secp256k1_pedersen_commitment_load(&tmpG, &coms[i]);
        secp256k1_gej_add_ge(&Cj, &hnegj, &tmpG);
        secp256k1_ge_set_gej(&tmpG, &Cj);
        secp256k1_pedersen_commitment_save(&C[i], &tmpG);
    }

    return secp256k1_unlinked_logarithmic_zero_com_verify(ctx, proof->data, ring_size, n, challenge, C, value_gen, blind_gen);
}



/*
 * const secp256k1_context* ctx,
        unsigned char *proof,
        int index,
        const unsigned char *blind,
        int ring_size,
        int n,
        const unsigned char *challenge,
        const secp256k1_pedersen_commitment *coms,
        const secp256k1_generator *value_gen,
        const secp256k1_generator *blind_gen
 */

int secp256k1_create_zero_mcom_short(const secp256k1_context* ctx,
                                     const ringcip_context *rctx,
                                     uint8_t *proof,
                                     int index,
                                     const unsigned char *blind,
                                     unsigned char *challenge,
                                     const secp256k1_pedersen_commitment *coms,
                                     int ring_size) {
    ARG_CHECK(ctx != NULL);
    ARG_CHECK(rctx != NULL);
    ARG_CHECK(coms != NULL);
    ARG_CHECK(blind != NULL);
    ARG_CHECK(index < ring_size);
    ARG_CHECK(ring_size <= rctx->N);

    uint8_t buf[32];
    secp256k1_scalar key;
    secp256k1_scalar rA, rB, rC, rD, rho[rctx->m], a[rctx->m][rctx->n], p[ring_size][rctx->m+1], tmpD, tmpA, tmpP;
    secp256k1_scalar x, f[rctx->m][rctx->n], z, zA, zC, tmp, tmp1;
    secp256k1_ge A, B, C, D, Q[rctx->m], tmpG, tmpH, tmpMu, tmpC;
    secp256k1_gej tmpj;
    int overflow, t, i, j;
    secp256k1_sha256 sha;

    secp256k1_generator_load(&tmpG, &rctx->geng);
    secp256k1_generator_load(&tmpH, &rctx->genh);
    secp256k1_generator_load(&tmpMu, &rctx->genmu);

    secp256k1_scalar_set_b32(&key, blind, &overflow);
    if (overflow) {
        return 0;
    }
    secp256k1_ecmult_const(&tmpj, &tmpG, &key, 256);
    secp256k1_ge_set_gej(&tmpC, &tmpj);
    uint8_t buf2[33];
    secp256k1_ge_save(buf2, &tmpC);
    if (memcmp(buf2, coms[index].data, 33) != 0) {
        printf("invalid key or commitment\n");
        return 0;
    }

    secp256k1_rand256(buf);
    secp256k1_scalar_set_b32(&rA, buf, &overflow);
    secp256k1_rand256(buf);
    secp256k1_scalar_set_b32(&rB, buf, &overflow);
    secp256k1_rand256(buf);
    secp256k1_scalar_set_b32(&rC, buf, &overflow);
    secp256k1_rand256(buf);
    secp256k1_scalar_set_b32(&rD, buf, &overflow);
    for (j = 0; j < rctx->m; j++) {
        secp256k1_rand256(buf);
        secp256k1_scalar_set_b32(&rho[j], buf, &overflow);
        secp256k1_scalar_set_int(&tmp, 0);
        for (i = 1; i < rctx->n; i++) {
            secp256k1_rand256(buf);
            secp256k1_scalar_set_b32(&a[j][i], buf, &overflow);
            secp256k1_scalar_add(&tmp, &tmp, &a[j][i]);
        }
        secp256k1_scalar_negate(&a[j][0], &tmp); // a[j][0] = - sum_{i=1}^{n} a[j][i]
    }
    secp256k1_scalar_set_b32(&key, blind, &overflow);
    if (overflow) {
        return 0;
    }

    // p_{t}(x) = prod_{j=0}^m (a_{j,t_j} + delta_{hat{t}_j, t_j}x)
    for (t = 0; t < ring_size; t++) {
        /* Set p[i][l] = 0 */
        for (j = 0; j < rctx->m + 1; j++) {
            secp256k1_scalar_set_int(&p[t][j], 0);
        }

        // copy first vector
        j = 0;
        //printf("%d %d %d ===\n",t, get_jth(j, rctx->n, t), get_jth(j, rctx->n, index));
        secp256k1_scalar_add(&p[t][j], &p[t][j], &a[j][get_jth(j, rctx->n, t)]);
        secp256k1_scalar_set_int(&p[t][j+1], get_jth(j, rctx->n, index) == get_jth(j, rctx->n, t));
        //start multiply
        for (j = 1; j < rctx->m; j++) {
            // set multiplier
            secp256k1_scalar_set_int(&tmpA, 0);
            secp256k1_scalar_add(&tmpA, &tmpA, &a[j][get_jth(j, rctx->n, t)]);
            secp256k1_scalar_set_int(&tmpD, get_jth(j, rctx->n, index) == get_jth(j, rctx->n, t));

            secp256k1_scalar_set_int(&tmpP, 0);
            secp256k1_scalar_add(&tmpP, &tmpP, &p[t][0]); // save original p[t][d]
            secp256k1_scalar_mul(&p[t][0], &p[t][0],&tmpA);
            for (int d = 1; d <= j; d++) {
                //printf("%d %d %d\n",t, j, d);
                secp256k1_scalar_mul(&tmp, &tmpP,&tmpD); // delta_{hat{t}_j, t_j} * p[t][d-1]

                secp256k1_scalar_set_int(&tmpP, 0);
                secp256k1_scalar_add(&tmpP, &tmpP, &p[t][d]); // save original p[t][d]

                secp256k1_scalar_mul(&p[t][d], &p[t][d],&tmpA); // a_{j,t_j} * p[t][d]
                secp256k1_scalar_add(&p[t][d], &p[t][d], &tmp); // p[t][d] = a_{j,t_j} * p[t][d] + delta_{hat{t}_j, t_j} * p[t][d-1]
            }
            secp256k1_scalar_mul(&p[t][j+1], &tmpP,&tmpD); // delta_{hat{t}_j, t_j} * p[t][d-1]
        }
        secp256k1_scalar_get_b32(buf, &p[t][rctx->m]);
    }

    // B := h^{r_B} prod_{j=0}^m prod_{i=0}^n h_{j,i}^{delta_{hat{t}_j,i}}
    secp256k1_ecmult_const(&tmpj, &tmpH, &rB, 256);
    secp256k1_ge_set_gej(&B, &tmpj);
    for (j = 0; j < rctx->m; j++)  {
        for (i = 0; i < rctx->n; i++) {
            secp256k1_scalar_set_int(&tmp, get_jth(j, rctx->n, index) == i);
            secp256k1_ecmult_const(&tmpj, &rctx->multigen[j*rctx->n + i], &tmp, 256);
            secp256k1_gej_add_ge(&tmpj, &tmpj, &B);
            secp256k1_ge_set_gej(&B, &tmpj);
        }
    }

    // A := h^{r_A} prod_{j=0}^m prod_{i=0}^n h_{j,i}^{a_{j,i}}
    secp256k1_ecmult_const(&tmpj, &tmpH, &rA, 256);
    secp256k1_ge_set_gej(&A, &tmpj);
    for (j = 0; j < rctx->m; j++)  {
        for (i = 0; i < rctx->n; i++) {
            secp256k1_ecmult_const(&tmpj, &rctx->multigen[j*rctx->n + i], &a[j][i], 256);
            secp256k1_gej_add_ge(&tmpj, &tmpj, &A);
            secp256k1_ge_set_gej(&A, &tmpj);
        }
    }

    // C := h^{r_C} prod_{j=0}^m prod_{i=0}^n h_{j,i}^{a_{j,i}(1-2 delta_{hat{t}_j,i})}
    secp256k1_ecmult_const(&tmpj, &tmpH, &rC, 256);
    secp256k1_ge_set_gej(&C, &tmpj);
    for (j = 0; j < rctx->m; j++)  {
        for (i = 0; i < rctx->n; i++) {
            secp256k1_scalar_set_int(&tmp, 2*(get_jth(j, rctx->n, index) == i));
            secp256k1_scalar_negate(&tmp, &tmp);
            secp256k1_scalar_set_int(&tmp1, 1);
            secp256k1_scalar_add(&tmp, &tmp, &tmp1);
            secp256k1_scalar_mul(&tmp, &tmp, &a[j][i]);
            secp256k1_ecmult_const(&tmpj, &rctx->multigen[j*rctx->n + i], &tmp, 256);
            secp256k1_gej_add_ge(&tmpj, &tmpj, &C);
            secp256k1_ge_set_gej(&C, &tmpj);
        }
    }

    // D := h^{r_D} prod_{j=0}^m prod_{i=0}^n h_{j,i}^{a^2_{j,i}}
    secp256k1_ecmult_const(&tmpj, &tmpH, &rD, 256);
    secp256k1_ge_set_gej(&D, &tmpj);
    for (j = 0; j < rctx->m; j++)  {
        for (i = 0; i < rctx->n; i++) {
            secp256k1_scalar_mul(&tmp, &a[j][i], &a[j][i]);
            secp256k1_scalar_negate(&tmp, &tmp);
            secp256k1_ecmult_const(&tmpj, &rctx->multigen[j*rctx->n + i], &tmp, 256);
            secp256k1_gej_add_ge(&tmpj, &tmpj, &D);
            secp256k1_ge_set_gej(&D, &tmpj);
        }
    }

    // Q_j := g^{rho_{j}}prod_{t=0}^N C_{t}^{p_{t,j}}
    for (j = 0; j < rctx->m; j++) {
        secp256k1_ecmult_const(&tmpj, &tmpG, &rho[j], 256);
        secp256k1_ge_set_gej(&Q[j], &tmpj);
        for (t = 0; t < ring_size; t++) {
            secp256k1_pedersen_commitment_load(&tmpC, &coms[t]);
            secp256k1_ecmult_const(&tmpj, &tmpC, &p[t][j], 256);
            secp256k1_gej_add_ge(&tmpj, &tmpj, &Q[j]);
            secp256k1_ge_set_gej(&Q[j], &tmpj);
        }
    }

    // Update proof
    int pointer = 0;
    secp256k1_ge_save(proof, &A);
    pointer += 33;
    secp256k1_ge_save(proof + pointer, &B);
    pointer += 33;
    secp256k1_ge_save(proof + pointer, &C);
    pointer += 33;
    secp256k1_ge_save(proof + pointer, &D);
    pointer += 33;
    for (j = 0; j < rctx->m; j++) {
        secp256k1_ge_save(proof + pointer, &Q[j]);
        pointer += 33;
    }

    // get x
    secp256k1_sha256_initialize(&sha);
    for (t = 0; t < ring_size; t++) {
        secp256k1_sha256_write(&sha, coms[t].data, 33);
    }
    secp256k1_sha256_write(&sha, proof, pointer);
    secp256k1_sha256_write(&sha, challenge, 32);
    secp256k1_sha256_finalize(&sha, buf);
    secp256k1_scalar_set_b32(&x, buf, NULL);

    // create f[j][i] = a[j][i] + delta_{hat{t}_j, i}x
    for (j = 0; j < rctx->m; j++) {
        for (i = 0; i < rctx->n; i++) {
            secp256k1_scalar_set_int(&f[j][i], get_jth(j, rctx->n, index) == i);
            secp256k1_scalar_mul(&f[j][i], &f[j][i], &x);
            secp256k1_scalar_add(&f[j][i], &f[j][i], &a[j][i]);
        }
    }

    // z_{A} := r_Bx + r_A
    secp256k1_scalar_mul(&tmp, &x, &rB);
    secp256k1_scalar_add(&zA, &tmp, &rA);

    // z_{C} := r_Cx + r_D
    secp256k1_scalar_mul(&tmp, &x, &rC);
    secp256k1_scalar_add(&zC, &tmp, &rD);

    // z := kx^m - sum_{j=0}^{m} rho_jx^j
    secp256k1_scalar_set_int(&z, 0);
    secp256k1_scalar_set_int(&tmp1, 1);
    for (j = 0; j < rctx->m; j++) {
        secp256k1_scalar_mul(&tmp, &rho[j], &tmp1);
        secp256k1_scalar_add(&z, &z, &tmp);
        secp256k1_scalar_mul(&tmp1, &tmp1, &x);
    }
    secp256k1_scalar_negate(&z, &z);
    secp256k1_scalar_mul(&tmp, &key, &tmp1);
    secp256k1_scalar_add(&z, &z, &tmp);

    // Update proof
    for (j = 0; j < rctx->m; j++) {
        for (i = 1; i < rctx->n; i++) {
            secp256k1_scalar_get_b32(proof + pointer, &f[j][i]);
            pointer += 32;
        }
    }
    secp256k1_scalar_get_b32(proof + pointer, &zA);
    pointer += 32;
    secp256k1_scalar_get_b32(proof + pointer, &zC);
    pointer += 32;
    secp256k1_scalar_get_b32(proof + pointer, &z);

    return 1;
}

int secp256k1_verify_zero_mcom_short(
        const secp256k1_context* ctx,
        const ringcip_context *rctx,
        uint8_t *proof,
        unsigned char *challenge,
        const secp256k1_pedersen_commitment *coms,
        int ring_size
        ) {
    ARG_CHECK(ctx != NULL);
    ARG_CHECK(rctx != NULL);
    ARG_CHECK(coms != NULL);
    ARG_CHECK(rctx->N >= ring_size);

    uint8_t buf[32];
    secp256k1_scalar x, f[rctx->m][rctx->n], z, zA, zC, tmp, tmp1, tmp2;
    secp256k1_ge A, B, C, D, Q[rctx->m], tmpG, tmpH, tmpMu, tmpC;
    secp256k1_gej tmpj;
    int overflow, t, i, j;
    secp256k1_sha256 sha;

    secp256k1_generator_load(&tmpG, &rctx->geng);
    secp256k1_generator_load(&tmpH, &rctx->genh);
    secp256k1_generator_load(&tmpMu, &rctx->genmu);

    // Load the proof
    int pointer = 0;
    secp256k1_ge_load(proof, &A);
    pointer += 33;
    secp256k1_ge_load(proof + pointer, &B);
    pointer += 33;
    secp256k1_ge_load(proof + pointer, &C);
    pointer += 33;
    secp256k1_ge_load(proof + pointer, &D);
    pointer += 33;
    for (j = 0; j < rctx->m; j++) {
        secp256k1_ge_load(proof + pointer, &Q[j]);
        pointer += 33;
    }

    // create x
    secp256k1_sha256_initialize(&sha);
    for (t = 0; t < ring_size; t++) {
        secp256k1_sha256_write(&sha, coms[t].data, 33);
    }
    secp256k1_sha256_write(&sha, proof, pointer);
    secp256k1_sha256_write(&sha, challenge, 32);
    secp256k1_sha256_finalize(&sha, buf);
    secp256k1_scalar_set_b32(&x, buf, NULL);

    // Load the rest of the proof
    for (j = 0; j < rctx->m; j++) {
        secp256k1_scalar_set_int(&tmp, 0);
        for (i = 1; i < rctx->n; i++) {
            secp256k1_scalar_set_b32(&f[j][i], proof + pointer, &overflow);
            pointer += 32;
            if (overflow) {
                return 0;
            }
            secp256k1_scalar_add(&tmp, &tmp, &f[j][i]);
        }
        secp256k1_scalar_negate(&f[j][0], &tmp); // f[j][0] = x - sum_{i=1}^{n} f[j][i]
        secp256k1_scalar_add(&f[j][0], &f[j][0], &x);
    }
    secp256k1_scalar_set_b32(&zA, proof + pointer, &overflow);
    if (overflow) {
        return 0;
    }
    pointer += 32;
    secp256k1_scalar_set_b32(&zC, proof + pointer, &overflow);
    if (overflow) {
        return 0;
    }
    pointer += 32;
    secp256k1_scalar_set_b32(&z, proof + pointer, &overflow);
    if (overflow) {
        return 0;
    }

    secp256k1_ge LHS;
    secp256k1_ge RHS;

    // prod_{t=0}^N C_{t}^{prod_{j=0}^m f_{j,t_j}} prod_{j=0}^m Q_j^{-x^j} stackrel{?}{=} g^{z}
    secp256k1_scalar_set_int(&tmp1, 1);
    for (j = 0; j < rctx->m; j++) {
        secp256k1_scalar_negate(&tmp2, &tmp1);
        secp256k1_ecmult_const(&tmpj, &Q[j], &tmp2, 256);
        if(j == 0) {
            secp256k1_ge_set_gej(&LHS, &tmpj);
        } else {
            secp256k1_gej_add_ge(&tmpj, &tmpj, &LHS);
            secp256k1_ge_set_gej(&LHS, &tmpj);
        }
        secp256k1_scalar_mul(&tmp1, &tmp1, &x);
    }

    for (t = 0; t < ring_size; t++) {
        secp256k1_scalar_set_int(&tmp, 1);
        for (j = 0; j < rctx->m; j++) {
            secp256k1_scalar_mul(&tmp, &tmp, &f[j][get_jth(j, rctx->n, t)]);
        }
        secp256k1_pedersen_commitment_load(&tmpC, &coms[t]);
        secp256k1_ecmult_const(&tmpj, &tmpC, &tmp, 256);
        secp256k1_gej_add_ge(&tmpj, &tmpj, &LHS);
        secp256k1_ge_set_gej(&LHS, &tmpj);
    }

    secp256k1_ecmult_const(&tmpj, &tmpG, &z, 256);
    secp256k1_ge_set_gej(&RHS, &tmpj);

    uint8_t lbuf[33];
    uint8_t rbuf[33];
    secp256k1_ge_save(lbuf, &LHS);
    secp256k1_ge_save(rbuf, &RHS);
    if (memcmp(lbuf, rbuf, 33) != 0) {
        return 0;
    }


    // B^xA stackrel{?}{=} h^{z_A} prod_{j=0}^m prod_{i=0}^n h_{j,i}^{f_{j,i}}
    secp256k1_ecmult_const(&tmpj, &B, &x, 256);
    secp256k1_gej_add_ge(&tmpj, &tmpj, &A);
    secp256k1_ge_set_gej(&LHS, &tmpj);

    secp256k1_ecmult_const(&tmpj, &tmpH, &zA, 256);
    secp256k1_ge_set_gej(&RHS, &tmpj);
    for (j = 0; j < rctx->m; j++)  {
        for (i = 0; i < rctx->n; i++) {
            secp256k1_ecmult_const(&tmpj, &rctx->multigen[j*rctx->n + i], &f[j][i], 256);
            secp256k1_gej_add_ge(&tmpj, &tmpj, &RHS);
            secp256k1_ge_set_gej(&RHS, &tmpj);
        }
    }

    secp256k1_ge_save(lbuf, &LHS);
    secp256k1_ge_save(rbuf, &RHS);
    if (memcmp(lbuf, rbuf, 33) != 0) {
        return 0;
    }

    // C^xD stackrel{?}{=} h^{z_C} prod_{j=0}^m prod_{i=0}^n h_{j,i}^{f_{j,i}(x-f_{j,i})}
    secp256k1_ecmult_const(&tmpj, &C, &x, 256);
    secp256k1_gej_add_ge(&tmpj, &tmpj, &D);
    secp256k1_ge_set_gej(&LHS, &tmpj);

    secp256k1_ecmult_const(&tmpj, &tmpH, &zC, 256);
    secp256k1_ge_set_gej(&RHS, &tmpj);
    for (j = 0; j < rctx->m; j++)  {
        for (i = 0; i < rctx->n; i++) {
            secp256k1_scalar_negate(&tmp, &f[j][i]);
            secp256k1_scalar_add(&tmp, &tmp, &x);
            secp256k1_scalar_mul(&tmp, &tmp, &f[j][i]);
            secp256k1_ecmult_const(&tmpj, &rctx->multigen[j*rctx->n + i], &tmp, 256);
            secp256k1_gej_add_ge(&tmpj, &tmpj, &RHS);
            secp256k1_ge_set_gej(&RHS, &tmpj);
        }
    }

    secp256k1_ge_save(lbuf, &LHS);
    secp256k1_ge_save(rbuf, &RHS);
    if (memcmp(lbuf, rbuf, 33) != 0) {
        return 0;
    }

    return 1;
}



#endif /* RAHAS_SECP256K1_MAIN_IMPL_H */