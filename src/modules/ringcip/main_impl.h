/**********************************************************************
 * Copyright (c) 2024 Jayamine Alupotha                                 *
 * Distributed under the MIT software license, see the accompanying   *
 * file COPYING2 or http://www.opensource.org/licenses/mit-license.php.*
 **********************************************************************/

#ifndef SECP256K_NOPENI_MAIN_IMPL_H
#define SECP256K_NOPENI_MAIN_IMPL_H

#include "include/secp256k1_ringcip.h"
#include "modules/nopenena/zero_com_impl.h"

static void secp256k1_ge_save1(unsigned char *buf, secp256k1_ge* ge) {
    secp256k1_fe_normalize(&ge->x);
    secp256k1_fe_get_b32(&buf[1], &ge->x);
    buf[0] = 9 ^ secp256k1_fe_is_quad_var(&ge->y);
}

static void secp256k1_ge_load1(unsigned char *buf, secp256k1_ge* ge) {
    secp256k1_fe fe;
    secp256k1_fe_set_b32(&fe, &buf[1]);
    secp256k1_ge_set_xquad(ge, &fe);
    if (buf[0] & 1) {
        secp256k1_ge_neg(ge, ge);
    }
}

ringcip_context secp256k1_ringcip_context_create(const secp256k1_context* ctx,
                                                 int L, int n, int m, uint8_t *gen_seed,
                                                 secp256k1_generator *blind_gen) {
    ringcip_context rctx;
    secp256k1_sha256 sha;
    unsigned char buf[32];
    secp256k1_generator gen;
    int j;
    rctx.L = L;
    rctx.n = n;
    rctx.m = m;
    rctx.N = 1;
    for (j = 0; j < rctx.m; j++) {
        rctx.N *= n;
    }
    memcpy(buf, gen_seed, 32);
    if (blind_gen == NULL) {
        while (1) {
            secp256k1_sha256_initialize(&sha);
            secp256k1_sha256_write(&sha, buf, 32);
            secp256k1_sha256_finalize(&sha, buf);

            if (secp256k1_generator_generate(ctx, &rctx.geng, buf))
                break;
        }
    }
    else {
        unsigned char gen_buffer[64];
        secp256k1_generator_serialize(ctx, gen_buffer, blind_gen);
        int ret = secp256k1_generator_parse(ctx, &rctx.geng, gen_buffer);
        if (ret == 0) {
            printf("bad blind gen in rctx\n");
        }
    }

    while (1) {
        secp256k1_sha256_initialize(&sha);
        secp256k1_sha256_write(&sha, buf, 32);
        secp256k1_sha256_finalize(&sha, buf);

        if(secp256k1_generator_generate(ctx, &rctx.genh, buf))
            break;
    }
    while (1) {
        secp256k1_sha256_initialize(&sha);
        secp256k1_sha256_write(&sha, buf, 32);
        secp256k1_sha256_finalize(&sha, buf);

        if(secp256k1_generator_generate(ctx, &rctx.genmu, buf))
            break;
    }
    rctx.multigen = (secp256k1_ge*) malloc(rctx.n*rctx.m*sizeof(secp256k1_ge));
    for (j = 0; j < rctx.n*rctx.m; j++) {
        while (1) {
            secp256k1_sha256_initialize(&sha);
            secp256k1_sha256_write(&sha, buf, 32);
            secp256k1_sha256_finalize(&sha, buf);

            if(secp256k1_generator_generate(ctx, &gen, buf))
                break;
        }
        secp256k1_generator_load(&rctx.multigen[j], &gen);
    }
    return rctx;
}

void secp256k1_ringcip_context_clear(ringcip_context *rctx) {
    free(rctx->multigen);
}

int secp256k1_create_secret_cint(const secp256k1_context* ctx, ringcip_context *rctx, cint_st *csk, int64_t v, uint8_t *key) {
    ARG_CHECK(ctx != NULL);
    ARG_CHECK(rctx != NULL);
    ARG_CHECK(csk != NULL);
    ARG_CHECK(key != NULL);
    int overflow;

    secp256k1_scalar_set_b32(&csk->key, key, &overflow);
    if (overflow) {
        return 0;
    }
    csk->v = v;

    return 1;
}

int secp256k1_create_cint(const secp256k1_context* ctx, ringcip_context *rctx, cint_pt *c, cint_st *csk) {
    ARG_CHECK(ctx != NULL);
    ARG_CHECK(rctx != NULL);
    ARG_CHECK(c != NULL);
    ARG_CHECK(csk != NULL);

    secp256k1_ge tmpG;
    secp256k1_ge tmpMu;
    secp256k1_gej tmpGj;

    /* g, mu */
    secp256k1_generator_load(&tmpG, &rctx->geng);
    secp256k1_generator_load(&tmpMu, &rctx->genmu);

    secp256k1_pedersen_ecmult(&tmpGj, &csk->key, csk->v, &tmpMu, &tmpG);
    secp256k1_ge_set_gej(&c->c, &tmpGj);

    // serialize
    secp256k1_ge_save1(c->buf, &c->c);

    return 1;
}


int secp256k1_serialize_cint(const secp256k1_context* ctx, uint8_t *buf, cint_pt *c) {
    ARG_CHECK(c != NULL);
    ARG_CHECK(buf!= NULL);
    memcpy(&c->buf, buf, 33);
    return 1;
}

int secp256k1_parse_cint(const secp256k1_context* ctx, cint_pt *c, uint8_t *buf) {
    ARG_CHECK(c != NULL);
    ARG_CHECK(buf!= NULL);
    secp256k1_ge_load1(buf, &c->c);
    memcpy(&c->buf, buf, 33);
    return 1;
}

int secp256k1_create_multival_com(const secp256k1_context* ctx, const ringcip_context *rctx, secp256k1_ge *com, secp256k1_scalar *vals, secp256k1_scalar *key) {
    ARG_CHECK(ctx != NULL);
    ARG_CHECK(rctx != NULL);
    ARG_CHECK(com != NULL);
    ARG_CHECK(key != NULL);
    int i;

    secp256k1_ge tmp;
    secp256k1_gej tmpj;

    /* g, mu */
    secp256k1_generator_load(&tmp, &rctx->genh);
    secp256k1_ecmult_const(&tmpj, &tmp, key, 256);
    secp256k1_ge_set_gej(com, &tmpj);
    for (i = 0; i < rctx->n * rctx->m; i++) {
        secp256k1_ecmult_const(&tmpj, &rctx->multigen[i], &vals[i], 256); // todo
        secp256k1_gej_add_ge(&tmpj, &tmpj, com);
        secp256k1_ge_set_gej(com, &tmpj);
    }
    return 1;
}

int get_jth(int j, int n, int num) {
    int i = 0, rem = num;
    while (i < j) {
        rem = rem / n;
        i++;
    }
    return rem - (rem/n)*n;
}


int secp256k1_zero_mcom_get_size(const ringcip_context *rctx) {
    return 33*(4+rctx->m) + 32*(3+(rctx->m*rctx->n));
}


int secp256k1_create_zero_mcom_proof(const secp256k1_context* ctx, const ringcip_context *rctx,
                                     uint8_t *proof, cint_pt *Cs, int index, secp256k1_scalar *key) {
    ARG_CHECK(ctx != NULL);
    ARG_CHECK(rctx != NULL);
    ARG_CHECK(Cs != NULL);
    ARG_CHECK(key != NULL);
    ARG_CHECK(index < rctx->N);

    uint8_t buf[32];
    secp256k1_scalar rA, rB, rC, rD, rho[rctx->m], a[rctx->m][rctx->n], p[rctx->N][rctx->m+1], tmpD, tmpA, tmpP;
    secp256k1_scalar x, f[rctx->m][rctx->n], z, zA, zC, tmp, tmp1;
    secp256k1_ge A, B, C, D, Q[rctx->m], tmpG, tmpH, tmpMu;
    secp256k1_gej tmpj;
    int overflow, t, i, j;
    secp256k1_sha256 sha;

    secp256k1_generator_load(&tmpG, &rctx->geng);
    secp256k1_generator_load(&tmpH, &rctx->genh);
    secp256k1_generator_load(&tmpMu, &rctx->genmu);

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
    if (overflow) {
        return 0;
    }

    // p_{t}(x) = prod_{j=0}^m (a_{j,t_j} + delta_{hat{t}_j, t_j}x)
    for (t = 0; t < rctx->N; t++) {
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
        for (t = 0; t < rctx->N; t++) {
            secp256k1_ecmult_const(&tmpj, &Cs[t].c, &p[t][j], 256);
            secp256k1_gej_add_ge(&tmpj, &tmpj, &Q[j]);
            secp256k1_ge_set_gej(&Q[j], &tmpj);
        }
    }

    // Update proof
    int pointer = 0;
    secp256k1_ge_save1(proof, &A);
    pointer += 33;
    secp256k1_ge_save1(proof + pointer, &B);
    pointer += 33;
    secp256k1_ge_save1(proof + pointer, &C);
    pointer += 33;
    secp256k1_ge_save1(proof + pointer, &D);
    pointer += 33;
    for (j = 0; j < rctx->m; j++) {
        secp256k1_ge_save1(proof + pointer, &Q[j]);
        pointer += 33;
    }

    // get x
    secp256k1_sha256_initialize(&sha);
    for (t = 0; t < rctx->N; t++) {
        secp256k1_sha256_write(&sha, Cs[t].buf, 33);
    }
    secp256k1_sha256_write(&sha, proof, pointer);
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
    secp256k1_scalar_mul(&tmp, key, &tmp1);
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


int secp256k1_verify_zero_mcom_proof(const secp256k1_context* ctx, const ringcip_context *rctx,
                                     uint8_t *proof, cint_pt *Cs) {
    ARG_CHECK(ctx != NULL);
    ARG_CHECK(rctx != NULL);
    ARG_CHECK(Cs != NULL);

    uint8_t buf[32];
    //secp256k1_scalar rA, rB, rC, rD, rho[rctx->m], a[rctx->m][rctx->n], p[rctx->N][rctx->m+1], tmpD, tmpA, tmpP;
    secp256k1_scalar x, f[rctx->m][rctx->n], z, zA, zC, tmp, tmp1, tmp2;
    secp256k1_ge A, B, C, D, Q[rctx->m], tmpG, tmpH, tmpMu;
    secp256k1_gej tmpj;
    int overflow, t, i, j;
    secp256k1_sha256 sha;

    secp256k1_generator_load(&tmpG, &rctx->geng);
    secp256k1_generator_load(&tmpH, &rctx->genh);
    secp256k1_generator_load(&tmpMu, &rctx->genmu);

    // Load the proof
    int pointer = 0;
    secp256k1_ge_load1(proof, &A);
    pointer += 33;
    secp256k1_ge_load1(proof + pointer, &B);
    pointer += 33;
    secp256k1_ge_load1(proof + pointer, &C);
    pointer += 33;
    secp256k1_ge_load1(proof + pointer, &D);
    pointer += 33;
    for (j = 0; j < rctx->m; j++) {
        secp256k1_ge_load1(proof + pointer, &Q[j]);
        pointer += 33;
    }

    // create x
    secp256k1_sha256_initialize(&sha);
    for (t = 0; t < rctx->N; t++) {
        secp256k1_sha256_write(&sha, Cs[t].buf, 33);
    }
    secp256k1_sha256_write(&sha, proof, pointer);
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

    for (t = 0; t < rctx->N; t++) {
        secp256k1_scalar_set_int(&tmp, 1);
        for (j = 0; j < rctx->m; j++) {
            secp256k1_scalar_mul(&tmp, &tmp, &f[j][get_jth(j, rctx->n, t)]);
        }
        secp256k1_ecmult_const(&tmpj, &Cs[t].c, &tmp, 256);
        secp256k1_gej_add_ge(&tmpj, &tmpj, &LHS);
        secp256k1_ge_set_gej(&LHS, &tmpj);
    }

    secp256k1_ecmult_const(&tmpj, &tmpG, &z, 256);
    secp256k1_ge_set_gej(&RHS, &tmpj);

    uint8_t lbuf[33];
    uint8_t rbuf[33];
    secp256k1_ge_save1(lbuf, &LHS);
    secp256k1_ge_save1(rbuf, &RHS);
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

    secp256k1_ge_save1(lbuf, &LHS);
    secp256k1_ge_save1(rbuf, &RHS);
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

    secp256k1_ge_save1(lbuf, &LHS);
    secp256k1_ge_save1(rbuf, &RHS);
    if (memcmp(lbuf, rbuf, 33) != 0) {
        return 0;
    }

    return 1;
}



#endif //SECP256K_NOPENI_MAIN_IMPL_H
