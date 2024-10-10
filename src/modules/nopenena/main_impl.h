/**********************************************************************
 * Copyright (c) 2024 Jayamine Alupotha                               *
 * Distributed under the MIT software license, see the accompanying   *
 * file COPYING2 or http://www.opensource.org/licenses/mit-license.php.*
 **********************************************************************/

#ifndef _SECP256K1_MODULE_NOPENI_
#define _SECP256K1_MODULE_NOPENI_

#include "include/secp256k1_nopenena.h"
#include "hash.h"
#include "math.h"
#include "zero_com_impl.h"



void scalar_print(char *msg, unsigned char *pbuf) {
    printf("%s:", msg);
    for (int i= 0; i < 10; i++) {
        printf("%d ", pbuf[i]);
    }
    printf("\n");
}

int get_ceiled_log_2(int N) {
    int t = 20; // only 2^12 max coms are supported
    int bit;
    int string;
    while (t >= 0) {
        string = 1 << t;
        bit = string & N;
        if (bit) {
            break;
        }
        t--;
    }
    if ((1 << t) < N) {
        t = t+1;
    }
    return t;
}

int factorial(int start, int end) {
    int fact = 1;
    while (start <= end) {
        fact *= start;
        start++;
    }
    return fact;
}

int compare(int *arr1, int *arr2, int p_size) {
    for (int i = 0; i < p_size; i++) {
        if (arr1[i] != arr2[i]) {
            return 0;
        }
    }
    return 1;
}

void get_combinations(int arr[], int data[], int *current, int *combinations, int start, int end, int index, int r) {
    if (index == r) {
        for (int j=0; j<r; j++) {
            combinations[(*current) * r + j] = data[j] - 1;
        }
        *current = *current + 1;
        return;
    }
    for (int i=start; i<=end && end-i+1 >= r-index; i++) {
        data[index] = arr[i];
        get_combinations(arr, data, current, combinations, i+1, end, index+1, r);
    }
}

int secp256k1_nopenena_get_j_of_all (int *j_indexes, int *combinations, int p_size, int total_combo) {
    for (int i = 0; i < total_combo; i++) {
        if (compare(combinations + (p_size*i), j_indexes, p_size)) {
            return i;
        }
    }
    return -1;
}

int secp256k1_nopenena_get_total_combo(int p_size, int n_size) {
    return factorial(n_size - p_size + 1, n_size)/ factorial(1, p_size);
}

int secp256k1_nopenena_get_allcombinations(int *combinations, int p_size, int n_size, int total_combo) {
    int arr[n_size];
    for (int i = 0; i < n_size; i++) {
        arr[i] = i + 1;
    }
    int data[p_size];
    int current_index = 0;
    get_combinations(arr, data, &current_index, combinations, 0, n_size-1, 0, p_size);
    return total_combo == current_index;
}

void secp256k1_open_nopenena_outputs_challenge(unsigned char *buf, nopenena_openout *out) {
    unsigned char gen_buf[33];
    secp256k1_sha256 sha;

    secp256k1_sha256_initialize(&sha);
    secp256k1_ge_save(gen_buf, &out->A_prime);
    secp256k1_sha256_write(&sha, gen_buf, 33);
    secp256k1_ge_save(gen_buf, &out->D);
    secp256k1_sha256_write(&sha, gen_buf, 33);
    for (int i = 0; i < out->p_size; i++) {
        secp256k1_ge_save(gen_buf, &out->C[i]);
        secp256k1_sha256_write(&sha, gen_buf, 33);
        secp256k1_ge_save(gen_buf, &out->A[i]);
        secp256k1_sha256_write(&sha, gen_buf, 33);
        secp256k1_ge_save(gen_buf, &out->B[i]);
        secp256k1_sha256_write(&sha, gen_buf, 33);

        secp256k1_scalar_get_b32(gen_buf, &out->f[i]);
        secp256k1_sha256_write(&sha, gen_buf, 32);
        secp256k1_scalar_get_b32(gen_buf, &out->z[i]);
        secp256k1_sha256_write(&sha, gen_buf, 32);
    }
    secp256k1_sha256_finalize(&sha, buf);
}

int native_log(int M, int n) {
    int i = 0;
    int M_ = 1;
    while (M > M_) {
        M_ *= n;
        i++;
    }
    return i;
}

void secp256k1_nopenena_context_create(
        const secp256k1_context* ctx,
        nopenena_context *nctx,
        int n, // for short zerocom
        int maximum_p_size,
        int maximum_acc_size,
        const unsigned char *gen_seed) {
    secp256k1_sha256 sha;
    unsigned char buf[32];
    int i;

    memset(buf, 0, 32);
    i = 0;
    while (i == 0) {
        secp256k1_sha256_initialize(&sha);
        secp256k1_sha256_write(&sha, gen_seed, 32);
        secp256k1_sha256_write(&sha, buf, 32);
        secp256k1_sha256_finalize(&sha, buf);

        if(secp256k1_generator_generate(ctx, &nctx->g, buf))
            break;
    }
    while (i == 0) {
        secp256k1_sha256_initialize(&sha);
        secp256k1_sha256_write(&sha, gen_seed, 32);
        secp256k1_sha256_write(&sha, buf, 32);
        secp256k1_sha256_finalize(&sha, buf);

        if(secp256k1_generator_generate(ctx, &nctx->mu, buf))
            break;
    }
    while (i == 0) {
        secp256k1_sha256_initialize(&sha);
        secp256k1_sha256_write(&sha, gen_seed, 32);
        secp256k1_sha256_write(&sha, buf, 32);
        secp256k1_sha256_finalize(&sha, buf);

        if(secp256k1_generator_generate(ctx, &nctx->h, buf))
            break;
    }
    unsigned char gen_seed1[32];
    memset(gen_seed1, 0, 32);
    int m = native_log(secp256k1_nopenena_get_total_combo(maximum_p_size, maximum_acc_size), n);
    //printf("vals: %d %d %d\n", m, maximum_p_size, maximum_acc_size);
    nctx->rctx = secp256k1_ringcip_context_create(ctx, 64, n, m, gen_seed1, &nctx->h);
    nctx->maximum_p_size = maximum_p_size;
    nctx->none = secp256k1_context_create(SECP256K1_CONTEXT_NONE);
    nctx->bp_gens = secp256k1_bulletproof_generators_create(nctx->none, &nctx->h, (maximum_p_size + 2)*2*64);;
    nctx->scratch = secp256k1_scratch_space_create(nctx->none, 1024 * 1024);
}

void secp256k1_nopenena_context_clear(const secp256k1_context* ctx, nopenena_context *nctx) {
    secp256k1_bulletproof_generators_destroy(ctx, nctx->bp_gens);
    secp256k1_scratch_destroy(nctx->scratch);
    secp256k1_context_destroy(nctx->none);
}

int secp256k1_create_account_sk(
        const secp256k1_context* ctx,
        nopenena_account_sk *acc_sk,
        unsigned char *r,
        unsigned char *k,
        unsigned char *gamma) {
    int overflow;
    ARG_CHECK(acc_sk != NULL);
    ARG_CHECK(r != NULL);
    ARG_CHECK(k != NULL);
    ARG_CHECK(gamma != NULL);

    secp256k1_scalar_set_b32(&acc_sk->r, r, &overflow);
    if (overflow) {
        return 0;
    }

    secp256k1_scalar_set_b32(&acc_sk->k, k, &overflow);
    if (overflow) {
        return 0;
    }

    secp256k1_scalar_set_b32(&acc_sk->gamma, gamma, &overflow);
    return overflow == 0;
}

int secp256k1_create_account_sk_r(
        const secp256k1_context* ctx,
        nopenena_account_sk *acc_sk,
        unsigned char *r) {
    int overflow;
    ARG_CHECK(acc_sk != NULL);
    ARG_CHECK(r != NULL);

    secp256k1_scalar_set_b32(&acc_sk->r, r, &overflow);
    return overflow == 0;
}

int secp256k1_create_openout_sk(
        const secp256k1_context* ctx,
        nopenena_openout_sk *openout_sk,
        unsigned char *alpha,
        unsigned char *a,
        unsigned char *a_prime,
        unsigned char *b,
        unsigned char *kappa,
        unsigned char *rho) {
    int overflow;
    ARG_CHECK(openout_sk != NULL);
    ARG_CHECK(alpha != NULL);
    ARG_CHECK(a != NULL);
    ARG_CHECK(a_prime != NULL);
    ARG_CHECK(b != NULL);
    ARG_CHECK(kappa != NULL);
    secp256k1_scalar_set_b32(&openout_sk->alpha, alpha, &overflow);
    if (overflow) {
        return 0;
    }
    secp256k1_scalar_set_b32(&openout_sk->a, a, &overflow);
    if (overflow) {
        return 0;
    }
    secp256k1_scalar_set_b32(&openout_sk->a_prime, a_prime, &overflow);
    if (overflow) {
        return 0;
    }
    secp256k1_scalar_set_b32(&openout_sk->b, b, &overflow);
    if (overflow) {
        return 0;
    }
    secp256k1_scalar_set_b32(&openout_sk->kappa, kappa, &overflow);
    if (overflow) {
        return 0;
    }

    secp256k1_scalar_set_b32(&openout_sk->rho, rho, &overflow);
    return overflow == 0;
}

int secp256k1_create_nopenena_account(
        const secp256k1_context* ctx,
        nopenena_context *nctx,
        nopenena_account *acc,
        nopenena_account_asset *asset,
        uint64_t value,
        nopenena_account_sk *acc_sk) {

    ARG_CHECK(acc != NULL);
    ARG_CHECK(asset != NULL);
    ARG_CHECK(acc_sk != NULL);
    ARG_CHECK(nctx != NULL);

    secp256k1_ge tmpG;
    secp256k1_ge tmpMu;
    secp256k1_gej tmpGj;

    /* g, mu */
    secp256k1_generator_load(&tmpG, &nctx->g);
    secp256k1_generator_load(&tmpMu, &nctx->mu);

    /* R = g^gamma */
    secp256k1_ecmult_const(&tmpGj, &tmpG, &acc_sk->gamma, 256);
    secp256k1_ge_set_gej(&acc->R, &tmpGj);

    /* K = g^{k*gamma} */
    secp256k1_ecmult_const(&tmpGj, &acc->R, &acc_sk->k, 256);
    secp256k1_ge_set_gej(&acc->K, &tmpGj);

    /* G = R^r */
    secp256k1_ecmult_const(&tmpGj, &acc->R, &acc_sk->r, 256);
    secp256k1_ge_set_gej(&asset->G, &tmpGj);

    /* V = K^r mu^v */
    secp256k1_pedersen_ecmult(&tmpGj, &acc_sk->r, value, &tmpMu, &acc->K);
    secp256k1_ge_set_gej(&asset->V, &tmpGj);

    return 1;
}


int secp256k1_update_nopenena_account(
        const secp256k1_context* ctx,
        nopenena_context *nctx,
        nopenena_account_asset *asset1,
        nopenena_account_asset *asset0,
        nopenena_account *acc,
        uint64_t value_new,
        uint64_t value,
        nopenena_account_sk *acc_sk) {

    ARG_CHECK(acc != NULL);
    ARG_CHECK(acc_sk != NULL);
    ARG_CHECK(asset1 != NULL);
    ARG_CHECK(asset0 != NULL);
    ARG_CHECK(nctx != NULL);

    secp256k1_ge tmpMu;
    secp256k1_gej tmpj;
    secp256k1_scalar tmp;
    secp256k1_scalar tmp1;

    /* mu */
    secp256k1_generator_load(&tmpMu, &nctx->mu);

    /* G = R^r */
    secp256k1_ecmult_const(&tmpj, &acc->R, &acc_sk->r, 256);
    secp256k1_gej_add_ge(&tmpj, &tmpj, &asset0->G);
    secp256k1_ge_set_gej(&asset1->G, &tmpj);

    /* V = K^r mu^v */
    secp256k1_ecmult_const(&tmpj, &acc->K, &acc_sk->r, 256);
    secp256k1_gej_add_ge(&tmpj, &tmpj, &asset0->V);
    secp256k1_ge_set_gej(&asset1->V, &tmpj);

    if ((value_new - value) != 0) {
        secp256k1_scalar_set_u64(&tmp, value);
        secp256k1_scalar_negate(&tmp, &tmp);
        secp256k1_scalar_set_u64(&tmp1, value_new);
        secp256k1_scalar_add(&tmp, &tmp, &tmp1);
        secp256k1_ecmult_const(&tmpj, &tmpMu, &tmp, 256);
        secp256k1_gej_add_ge(&tmpj, &tmpj, &asset1->V);
        secp256k1_ge_set_gej(&asset1->V, &tmpj);
    }

    return 1;
}



int secp256k1_prove_nopenena_update(
        const secp256k1_context* ctx,
        nopenena_context *nctx,
        nopenena_update_proof *proof,
        nopenena_account_asset *asset1,
        nopenena_account_asset *asset0,
        nopenena_account *acc,
        uint64_t value_new,
        uint64_t value,
        nopenena_account_sk *acc_sk,
        nopenena_openout_sk *openouts_sk,
        nopenena_openout *w,
        int flag,
        int i,
        const unsigned char *nonce) {
    ARG_CHECK(ctx != NULL);
    ARG_CHECK(nctx != NULL);
    ARG_CHECK(w != NULL);
    ARG_CHECK(acc != NULL);
    if (flag == 1) { // owner update
        ARG_CHECK(acc_sk != NULL);
        ARG_CHECK(openouts_sk != NULL);
        ARG_CHECK(asset1 != NULL);
        ARG_CHECK(asset0 != NULL);
        ARG_CHECK(value_new - value != 0);
    }
    secp256k1_scalar t;
    secp256k1_scalar tau;
    secp256k1_ge T1;
    secp256k1_ge T2;
    secp256k1_ge T3;
    secp256k1_sha256 sha;
    secp256k1_ge tmpMu;
    secp256k1_gej tmpj;
    secp256k1_gej tmpj1;
    unsigned char gen_buf[33];
    secp256k1_scalar x;
    secp256k1_scalar tmp;
    secp256k1_scalar tmp1;

    secp256k1_generator_load(&tmpMu, &nctx->mu);

    secp256k1_scalar_chacha20(&t, &tau, nonce, i);

    // T1 = R^t
    secp256k1_ecmult_const(&tmpj, &acc->R, &t, 256);
    secp256k1_ge_set_gej(&T1, &tmpj);
    // T2 = K^tmu^{tau}
    secp256k1_ecmult_const(&tmpj, &acc->K, &t, 256);
    secp256k1_ecmult_const(&tmpj1, &tmpMu, &tau, 256);
    secp256k1_ge_set_gej(&T2, &tmpj);
    secp256k1_gej_add_ge(&tmpj, &tmpj1, &T2);
    secp256k1_ge_set_gej(&T2, &tmpj);
    // T3 = K^{tau}R^{kappa}
    secp256k1_ecmult_const(&tmpj, &acc->K, &tau, 256);
    secp256k1_ecmult_const(&tmpj1, &acc->R, &openouts_sk->kappa, 256);
    secp256k1_ge_set_gej(&T3, &tmpj);
    secp256k1_gej_add_ge(&tmpj, &tmpj1, &T3);
    secp256k1_ge_set_gej(&T3, &tmpj);

    // x
    secp256k1_sha256_initialize(&sha);
    secp256k1_open_nopenena_outputs_challenge(gen_buf, w);
    secp256k1_sha256_write(&sha, gen_buf, 32); // 32 not 33
    secp256k1_ge_save(gen_buf, &acc->R);
    secp256k1_sha256_write(&sha, gen_buf, 33);
    secp256k1_ge_save(gen_buf, &acc->K);
    secp256k1_sha256_write(&sha, gen_buf, 33);
    secp256k1_ge_save(gen_buf, &asset0->G);
    secp256k1_sha256_write(&sha, gen_buf, 33);
    secp256k1_ge_save(gen_buf, &asset0->V);
    secp256k1_sha256_write(&sha, gen_buf, 33);
    secp256k1_ge_save(gen_buf, &asset1->G);
    secp256k1_sha256_write(&sha, gen_buf, 33);
    secp256k1_ge_save(gen_buf, &asset1->V);
    secp256k1_sha256_write(&sha, gen_buf, 33);
    secp256k1_ge_save(gen_buf, &T1);
    secp256k1_sha256_write(&sha, gen_buf, 33);
    secp256k1_ge_save(gen_buf, &T2);
    secp256k1_sha256_write(&sha, gen_buf, 33);
    secp256k1_ge_save(gen_buf, &T3);
    secp256k1_sha256_write(&sha, gen_buf, 33);
    secp256k1_sha256_finalize(&sha, proof->x);
    secp256k1_scalar_set_b32(&x, proof->x, NULL);

    // s_1
    secp256k1_scalar_mul(&tmp, &x, &acc_sk->r);
    secp256k1_scalar_add(&proof->s1, &tmp, &t);

    // s_2
    if (flag == 1) {
        secp256k1_scalar_set_u64(&tmp, value);
        secp256k1_scalar_negate(&tmp, &tmp);
        secp256k1_scalar_set_u64(&tmp1, value_new);
        secp256k1_scalar_add(&tmp, &tmp, &tmp1);
        secp256k1_scalar_mul(&tmp, &x, &tmp);
    } else {
        secp256k1_scalar_set_int(&tmp, 0);
    }
    secp256k1_scalar_add(&proof->s2, &tmp, &tau);

    // s_3
    if (flag == 1) {
        secp256k1_scalar_mul(&tmp, &acc_sk->k, &tmp);
    } else {
        secp256k1_scalar_set_int(&tmp, 0);
    }
    secp256k1_scalar_negate(&proof->s3, &openouts_sk->kappa);
    secp256k1_scalar_add(&proof->s3, &proof->s3, &tmp);

    return 1;

}

int secp256k1_verify_nopenena_update(
        const secp256k1_context* ctx,
        nopenena_context *nctx,
        nopenena_update_proof *proof,
        nopenena_account_asset *asset1,
        nopenena_account_asset *asset0,
        nopenena_account *acc,
        nopenena_openout *w) {
    ARG_CHECK(ctx != NULL);
    ARG_CHECK(nctx != NULL);
    ARG_CHECK(w != NULL);
    ARG_CHECK(acc != NULL);


    secp256k1_sha256 sha;
    secp256k1_ge tmpMu;
    secp256k1_ge tmpG;
    secp256k1_gej tmpj;
    unsigned char buf[32];
    unsigned char gen_buf[33];
    secp256k1_scalar x;
    secp256k1_scalar tmp;
    secp256k1_ge T10;
    secp256k1_ge T20;
    secp256k1_ge T30;

    secp256k1_generator_load(&tmpMu, &nctx->mu);
    secp256k1_scalar_set_b32(&x, proof->x, NULL);

    // T_1 = R^{s_1}((G')^{-1}G)^x
    secp256k1_ge_neg(&tmpG, &asset1->G);
    secp256k1_gej_set_ge(&tmpj, &tmpG);
    secp256k1_gej_add_ge(&tmpj, &tmpj, &asset0->G);
    secp256k1_ge_set_gej(&T10, &tmpj);
    secp256k1_ecmult_const(&tmpj, &T10, &x, 256);
    secp256k1_ge_set_gej(&T10, &tmpj);

    secp256k1_ecmult_const(&tmpj, &acc->R, &proof->s1, 256);
    secp256k1_gej_add_ge(&tmpj, &tmpj, &T10);
    secp256k1_ge_set_gej(&T10, &tmpj);

    // T_2 = K^{s_1}\mu^{s_2}((V')^{-1}V)^x
    secp256k1_ge_neg(&tmpG, &asset1->V);
    secp256k1_gej_set_ge(&tmpj, &tmpG);
    secp256k1_gej_add_ge(&tmpj, &tmpj, &asset0->V);
    secp256k1_ge_set_gej(&T20, &tmpj);
    secp256k1_ecmult_const(&tmpj, &T20, &x, 256);
    secp256k1_ge_set_gej(&T20, &tmpj);

    secp256k1_ecmult_const(&tmpj, &acc->K, &proof->s1, 256);
    secp256k1_gej_add_ge(&tmpj, &tmpj, &T20);
    secp256k1_ge_set_gej(&T20, &tmpj);

    secp256k1_ecmult_const(&tmpj, &tmpMu, &proof->s2, 256);
    secp256k1_gej_add_ge(&tmpj, &tmpj, &T20);
    secp256k1_ge_set_gej(&T20, &tmpj);


    // T_3 = K^{s_2}R^{-s_3}
    secp256k1_scalar_negate(&tmp, &proof->s3);
    secp256k1_ecmult_const(&tmpj, &acc->R, &tmp, 256);
    secp256k1_ge_set_gej(&T30, &tmpj);
    secp256k1_ecmult_const(&tmpj, &acc->K, &proof->s2, 256);
    secp256k1_gej_add_ge(&tmpj, &tmpj, &T30);
    secp256k1_ge_set_gej(&T30, &tmpj);


    // x
    secp256k1_sha256_initialize(&sha);
    secp256k1_open_nopenena_outputs_challenge(gen_buf, w);
    secp256k1_sha256_write(&sha, gen_buf, 32); // 32 not 33
    secp256k1_ge_save(gen_buf, &acc->R);
    secp256k1_sha256_write(&sha, gen_buf, 33);
    secp256k1_ge_save(gen_buf, &acc->K);
    secp256k1_sha256_write(&sha, gen_buf, 33);
    secp256k1_ge_save(gen_buf, &asset0->G);
    secp256k1_sha256_write(&sha, gen_buf, 33);
    secp256k1_ge_save(gen_buf, &asset0->V);
    secp256k1_sha256_write(&sha, gen_buf, 33);
    secp256k1_ge_save(gen_buf, &asset1->G);
    secp256k1_sha256_write(&sha, gen_buf, 33);
    secp256k1_ge_save(gen_buf, &asset1->V);
    secp256k1_sha256_write(&sha, gen_buf, 33);
    secp256k1_ge_save(gen_buf, &T10);
    secp256k1_sha256_write(&sha, gen_buf, 33);
    secp256k1_ge_save(gen_buf, &T20);
    secp256k1_sha256_write(&sha, gen_buf, 33);
    secp256k1_ge_save(gen_buf, &T30);
    secp256k1_sha256_write(&sha, gen_buf, 33);
    secp256k1_sha256_finalize(&sha, buf);

    return memcmp(buf, proof->x, 32) == 0;
}

int secp256k1_open_nopenena_outputs2(
        const secp256k1_context* ctx,
        nopenena_context *nctx,
        nopenena_openout *out,
        nopenena_account_asset *asset[],
        uint64_t *values_new,
        uint64_t *values,
        nopenena_account_sk *acc_sk[],
        nopenena_openout_sk *openouts_sk[],
        int *j_indexes,
        int p_size,
        int n_size,
        unsigned char *nonce) {

    ARG_CHECK(ctx != NULL);
    ARG_CHECK(nctx != NULL);
    ARG_CHECK(out != NULL);
    ARG_CHECK(asset != NULL);
    ARG_CHECK(values_new != NULL);
    ARG_CHECK(acc_sk != NULL);
    ARG_CHECK(openouts_sk != NULL);
    ARG_CHECK(nctx->maximum_p_size >= p_size);

    secp256k1_sha256 sha;
    secp256k1_ge tmpG;
    secp256k1_ge tmpMu;
    secp256k1_ge tmpH;
    secp256k1_gej tmpj[p_size];
    secp256k1_gej tmpDj;
    const unsigned char *blind_ptr[p_size];
    unsigned char buf[32];
    unsigned char gen_buf[33];
    unsigned char blinds[p_size][32];
    int i;
    secp256k1_scalar tmp;
    secp256k1_scalar tmp1;
    out->pi_len = sizeof(out->pi_range);

    out->p_size = p_size;
    /* g, h, mu */
    secp256k1_generator_load(&tmpG, &nctx->g);
    secp256k1_generator_load(&tmpH, &nctx->h);
    secp256k1_generator_load(&tmpMu, &nctx->mu);

    out->C = (secp256k1_ge *) malloc(p_size * sizeof(secp256k1_ge));
    out->C_bin = (secp256k1_pedersen_commitment *) malloc(p_size * sizeof(secp256k1_pedersen_commitment));

    // create assets
    for (i = 0; i < p_size; i++) {
        secp256k1_pedersen_ecmult(&tmpj[i], &openouts_sk[j_indexes[i]]->alpha, values_new[i], &tmpMu, &tmpH);
        secp256k1_ge_set_gej(&out->C[i], &tmpj[i]);
        secp256k1_pedersen_commitment_save(&out->C_bin[i], &out->C[i]);

        secp256k1_scalar_get_b32(blinds[i], &openouts_sk[j_indexes[i]]->alpha);
        blind_ptr[i] = blinds[i];
    }
    // range proofs
    int ret = secp256k1_bulletproof_rangeproof_prove(
            ctx,
            nctx->scratch,
            nctx->bp_gens,
            out->pi_range, &out->pi_len,
            NULL, NULL, NULL,
            values_new, NULL, blind_ptr, NULL, p_size, &nctx->mu, 64, nonce, NULL, NULL, 0, NULL);
    if (!ret) {
        free(out->C);
        free(out->C_bin);
        return 0;
    }

    out->A = (secp256k1_ge *) malloc(p_size * sizeof(secp256k1_ge));
    out->B = (secp256k1_ge *) malloc(p_size * sizeof(secp256k1_ge));
    out->f = (secp256k1_scalar *) malloc(p_size * sizeof(secp256k1_scalar));
    out->z = (secp256k1_scalar *) malloc(p_size * sizeof(secp256k1_scalar));
    out->y_a = (secp256k1_scalar *) malloc(p_size * sizeof(secp256k1_scalar));
    out->y_b = (secp256k1_scalar *) malloc(p_size * sizeof(secp256k1_scalar));

    for (i = 0; i < p_size; i++) {
        // A = C^{a}
        secp256k1_ecmult_const(&tmpj[i], &out->C[i], &openouts_sk[j_indexes[i]]->a, 256);
        secp256k1_ge_set_gej(&out->A[i], &tmpj[i]);

        // A' = h^{a'}
        secp256k1_ecmult_const(&tmpj[i], &tmpH, &openouts_sk[j_indexes[i]]->a_prime, 256);
        if (i==0) {
            secp256k1_ge_set_gej(&out->A_prime, &tmpj[i]);
        } else {
            secp256k1_gej_add_ge(&tmpj[i], &tmpj[i], &out->A_prime);
            secp256k1_ge_set_gej(&out->A_prime, &tmpj[i]);
        }

        // B = G'^{b}
        secp256k1_ecmult_const(&tmpj[i], &asset[j_indexes[i]]->G, &openouts_sk[j_indexes[i]]->b, 256);
        secp256k1_ge_set_gej(&out->B[i], &tmpj[i]);

        /* y_a */
        secp256k1_sha256_initialize(&sha);
        secp256k1_sha256_write(&sha, out->C_bin[i].data, 33);
        secp256k1_ge_save(gen_buf, &out->A[i]);
        secp256k1_sha256_write(&sha, gen_buf, 33);
        secp256k1_ge_save(gen_buf, &out->B[i]);
        secp256k1_sha256_write(&sha, gen_buf, 33);
        secp256k1_sha256_finalize(&sha, buf);
        secp256k1_scalar_set_b32(&out->y_a[i], buf, NULL);

        // fl := ρkl(v′l − vl) − ya,lal
        secp256k1_scalar_set_u64(&tmp, values[i]);
        secp256k1_scalar_negate(&tmp, &tmp);
        secp256k1_scalar_set_u64(&tmp1, values_new[i]);
        secp256k1_scalar_add(&tmp, &tmp, &tmp1);
        secp256k1_scalar_mul(&tmp, &tmp, &acc_sk[i]->k);
        secp256k1_scalar_mul(&tmp1, &tmp, &openouts_sk[j_indexes[i]]->rho);
        secp256k1_scalar_mul(&tmp, &out->y_a[i], &openouts_sk[j_indexes[i]]->a);
        secp256k1_scalar_negate(&tmp, &tmp);
        secp256k1_scalar_add(&out->f[i], &tmp, &tmp1);

        /* y_b */
        secp256k1_sha256_initialize(&sha);
        secp256k1_sha256_write(&sha, buf, 32);
        secp256k1_scalar_get_b32(buf, &out->f[i]);
        secp256k1_sha256_finalize(&sha, buf);
        secp256k1_scalar_set_b32(&out->y_b[i], buf, NULL);

        // zl
        secp256k1_scalar_mul(&tmp1, &tmp1, &acc_sk[i]->k);
        secp256k1_scalar_mul(&tmp, &out->y_b[i], &openouts_sk[j_indexes[i]]->b);
        secp256k1_scalar_negate(&tmp, &tmp);
        secp256k1_scalar_add(&out->z[i], &tmp1, &tmp);
    }

    // create D values
    secp256k1_ecmult_const(&tmpDj, &asset[0]->V, &openouts_sk[0]->kappa, 256);
    secp256k1_ge_set_gej(&out->D, &tmpDj);
    for (i = 1; i < n_size; i++) {
        secp256k1_ecmult_const(&tmpDj, &asset[i]->V, &openouts_sk[i]->kappa, 256);
        secp256k1_gej_add_ge(&tmpDj, &tmpDj, &out->D);
        secp256k1_ge_set_gej(&out->D, &tmpDj);
    }

    return 1;
}


int contract_escrow_pcompile(nopenena_context *nctx, nopenena_escrow_zk *f_zk,
                             unsigned char *alpha_e, unsigned char *alpha_cS, unsigned char *alpha_cB,
                             uint64_t e, uint64_t cS, uint64_t cB) {
    secp256k1_gej tmpj;
    secp256k1_scalar tmp;
    secp256k1_ge tmpMu;
    secp256k1_ge tmpH;
    secp256k1_ge C;

    secp256k1_generator_load(&tmpH, &nctx->h);
    secp256k1_generator_load(&tmpMu, &nctx->mu);

    secp256k1_scalar_set_b32(&tmp, alpha_e, NULL);
    secp256k1_pedersen_ecmult(&tmpj, &tmp, e, &tmpMu, &tmpH);
    secp256k1_ge_set_gej(&C, &tmpj);
    secp256k1_pedersen_commitment_save(&f_zk->E, &C);

    secp256k1_scalar_set_b32(&tmp, alpha_cS, NULL);
    secp256k1_pedersen_ecmult(&tmpj, &tmp, cS, &tmpMu, &tmpH);
    secp256k1_ge_set_gej(&C, &tmpj);
    secp256k1_pedersen_commitment_save(&f_zk->CS, &C);

    secp256k1_scalar_set_b32(&tmp, alpha_cB, NULL);
    secp256k1_pedersen_ecmult(&tmpj, &tmp, cB, &tmpMu, &tmpH);
    secp256k1_ge_set_gej(&C, &tmpj);
    secp256k1_pedersen_commitment_save(&f_zk->CB, &C);

    return 1;
}


int contract_escrow_prove(nopenena_context *nctx, secp256k1_pedersen_commitment *out_C_bin, unsigned char *pi_contract,
                          nopenena_escrow_zk *f_zk, unsigned char *alpha_e, unsigned char *alpha_c,
                          secp256k1_scalar *alpha_c0, secp256k1_scalar *alpha_c1,
                          unsigned char *y_1, unsigned char *y_2) {
    secp256k1_gej tmpj;
    secp256k1_scalar tmp;
    secp256k1_scalar s;
    secp256k1_scalar x;
    secp256k1_scalar y1;
    secp256k1_scalar y2;
    secp256k1_ge tmpMu;
    secp256k1_ge tmpH;
    secp256k1_ge R;
    unsigned char buf[32];

    int pointer = 0;

    secp256k1_generator_load(&tmpH, &nctx->h);
    secp256k1_generator_load(&tmpMu, &nctx->mu);

    // pi_eqE
    secp256k1_scalar_set_b32(&y1, y_1, NULL);
    secp256k1_ecmult_const(&tmpj, &tmpH, &y1, 256);
    secp256k1_ge_set_gej(&R, &tmpj);
    secp256k1_ge_save(pi_contract, &R);
    pointer+= 33;

    // x
    secp256k1_sha256 sha;
    secp256k1_sha256_initialize(&sha);
    secp256k1_sha256_write(&sha, f_zk->E.data, 33);
    secp256k1_sha256_write(&sha, out_C_bin[1].data, 33);
    secp256k1_sha256_write(&sha, pi_contract, 33);
    secp256k1_sha256_finalize(&sha, buf);
    secp256k1_scalar_set_b32(&x, buf, NULL);

    // s
    secp256k1_scalar_negate(&s, alpha_c1);
    secp256k1_scalar_set_b32(&tmp, alpha_e, NULL);
    secp256k1_scalar_add(&s, &s, &tmp);
    secp256k1_scalar_mul(&s, &s, &x);
    secp256k1_scalar_add(&s, &s, &y1);
    secp256k1_scalar_get_b32(pi_contract + pointer, &s);
    pointer += 32;

    secp256k1_scalar_set_b32(&y2, y_2, NULL);
    secp256k1_ecmult_const(&tmpj, &tmpH, &y2, 256);
    secp256k1_ge_set_gej(&R, &tmpj);
    secp256k1_ge_save(pi_contract + pointer, &R);
    pointer += 33;

    // x
    secp256k1_sha256_initialize(&sha);
    secp256k1_sha256_write(&sha, f_zk->E.data, 33);
    secp256k1_sha256_write(&sha, f_zk->CS.data, 33);
    secp256k1_sha256_write(&sha, f_zk->CB.data, 33);
    secp256k1_sha256_write(&sha, out_C_bin[0].data, 33);
    secp256k1_sha256_write(&sha, pi_contract + 65, 33);
    secp256k1_sha256_finalize(&sha, buf);
    secp256k1_scalar_set_b32(&x, buf, NULL);

    // s
    secp256k1_scalar_negate(&s, alpha_c0);
    secp256k1_scalar_set_b32(&tmp, alpha_c, NULL);
    secp256k1_scalar_add(&s, &s, &tmp);
    secp256k1_scalar_mul(&s, &s, &x);
    secp256k1_scalar_add(&s, &s, &y2);
    secp256k1_scalar_get_b32(pi_contract + pointer, &s);
    pointer += 32;

    return 1;
}

int contract_escrow_verify(nopenena_context *nctx, secp256k1_pedersen_commitment *out_C_bin, unsigned char *pi_contract,
                           nopenena_escrow_zk *f_zk, int refund) {

    secp256k1_gej tmpj;
    secp256k1_scalar s;
    secp256k1_scalar x;
    secp256k1_ge tmpMu;
    secp256k1_ge tmpH;
    secp256k1_ge C1;
    secp256k1_ge C2;
    secp256k1_ge R;
    secp256k1_ge LHS;
    secp256k1_ge RHS;
    unsigned char lhs[33];
    unsigned char rhs[33];
    unsigned char buf[32];

    int pointer = 0;

    secp256k1_generator_load(&tmpH, &nctx->h);
    secp256k1_generator_load(&tmpMu, &nctx->mu);

    secp256k1_ge_load(pi_contract + pointer, &R);
    pointer+= 33;
    secp256k1_scalar_set_b32(&s, pi_contract + pointer, NULL);
    pointer+= 32;

    // x
    secp256k1_sha256 sha;
    secp256k1_sha256_initialize(&sha);
    secp256k1_sha256_write(&sha, f_zk->E.data, 33);
    secp256k1_sha256_write(&sha, out_C_bin[1].data, 33);
    secp256k1_sha256_write(&sha, pi_contract, 33);
    secp256k1_sha256_finalize(&sha, buf);
    secp256k1_scalar_set_b32(&x, buf, NULL);

    // s
    secp256k1_ecmult_const(&tmpj, &tmpH, &s, 256);
    secp256k1_ge_set_gej(&LHS, &tmpj);
    secp256k1_ge_save(lhs, &LHS);

    secp256k1_pedersen_commitment_load(&C1, &f_zk->E);
    secp256k1_pedersen_commitment_load(&C2, &out_C_bin[1]);
    secp256k1_ge_neg(&C2, &C2);
    secp256k1_gej_set_ge(&tmpj, &C1);
    secp256k1_gej_add_ge(&tmpj, &tmpj, &C2);
    secp256k1_ge_set_gej(&RHS, &tmpj);
    secp256k1_ecmult_const(&tmpj, &RHS, &x, 256);
    secp256k1_gej_add_ge(&tmpj,&tmpj, &R);
    secp256k1_ge_set_gej(&RHS, &tmpj);
    secp256k1_ge_save(rhs, &RHS);

    if (memcmp(rhs, lhs, 33) != 0)  {
        return 0;
    }

    if (refund == 0) {
        secp256k1_pedersen_commitment_load(&C1, &f_zk->CS);
        secp256k1_pedersen_commitment_load(&C2, &out_C_bin[0]);
    }
    if (refund == 1) {
        secp256k1_pedersen_commitment_load(&C1, &f_zk->CB);
        secp256k1_pedersen_commitment_load(&C2, &out_C_bin[0]);
    }

    secp256k1_ge_load(pi_contract + pointer, &R);
    pointer+= 33;
    secp256k1_scalar_set_b32(&s, pi_contract + pointer, NULL);
    pointer+= 32;

    // x
    secp256k1_sha256_initialize(&sha);
    secp256k1_sha256_write(&sha, f_zk->E.data, 33);
    secp256k1_sha256_write(&sha, f_zk->CS.data, 33);
    secp256k1_sha256_write(&sha, f_zk->CB.data, 33);
    secp256k1_sha256_write(&sha, out_C_bin[0].data, 33);
    secp256k1_sha256_write(&sha, pi_contract + 65, 33);
    secp256k1_sha256_finalize(&sha, buf);
    secp256k1_scalar_set_b32(&x, buf, NULL);

    // s
    secp256k1_ecmult_const(&tmpj, &tmpH, &s, 256);
    secp256k1_ge_set_gej(&LHS, &tmpj);
    secp256k1_ge_save(lhs, &LHS);

    secp256k1_ge_neg(&C2, &C2);
    secp256k1_gej_set_ge(&tmpj, &C1);
    secp256k1_gej_add_ge(&tmpj, &tmpj, &C2);
    secp256k1_ge_set_gej(&RHS, &tmpj);
    secp256k1_ecmult_const(&tmpj, &RHS, &x, 256);
    secp256k1_gej_add_ge(&tmpj,&tmpj, &R);
    secp256k1_ge_set_gej(&RHS, &tmpj);
    secp256k1_ge_save(rhs, &RHS);


    return memcmp(rhs, lhs, 33) == 0;
}

int secp256k1_open_nopenena_outputs(
        const secp256k1_context* ctx,
        nopenena_context *nctx,
        nopenena_openout *out,
        nopenena_account_asset *asset[],
        uint64_t *values_new,
        uint64_t *values,
        nopenena_account_sk *acc_sk[],
        nopenena_openout_sk *openouts_sk[],
        int *j_indexes,
        int p_size,
        int n_size,
        unsigned char *nonce,
        uint64_t withheld_c,
        unsigned char *alpha_c,
        uint64_t withheld_c_new,
        unsigned char *alpha_c_new) {

    ARG_CHECK(ctx != NULL);
    ARG_CHECK(nctx != NULL);
    ARG_CHECK(out != NULL);
    ARG_CHECK(asset != NULL);
    ARG_CHECK(values_new != NULL);
    ARG_CHECK(acc_sk != NULL);
    ARG_CHECK(openouts_sk != NULL);
    ARG_CHECK(nctx->maximum_p_size >= p_size);

    secp256k1_sha256 sha;
    secp256k1_ge tmpG;
    secp256k1_ge tmpMu;
    secp256k1_ge tmpH;
    secp256k1_gej tmpj[p_size+2];
    secp256k1_gej tmpDj;
    unsigned char buf[32];
    unsigned char gen_buf[33];
    int i;
    secp256k1_scalar tmp;
    secp256k1_scalar tmp1;
    out->pi_len = sizeof(out->pi_range);

    out->p_size = p_size;
    /* g, h, mu */
    secp256k1_generator_load(&tmpG, &nctx->g);
    secp256k1_generator_load(&tmpH, &nctx->h);
    secp256k1_generator_load(&tmpMu, &nctx->mu);

    out->range_proof_size = p_size;
    if (alpha_c != NULL)
        out->range_proof_size++;
    if (alpha_c_new != NULL)
        out->range_proof_size++;
    unsigned char blinds[out->range_proof_size][32];
    const unsigned char *blind_ptr[out->range_proof_size];
    out->C = (secp256k1_ge *) malloc(out->range_proof_size * sizeof(secp256k1_ge));
    out->C_bin = (secp256k1_pedersen_commitment *) malloc(out->range_proof_size * sizeof(secp256k1_pedersen_commitment));
    uint64_t range_proof_values[out->range_proof_size];


    // create assets
    for (i = 0; i < p_size; i++) {
        secp256k1_pedersen_ecmult(&tmpj[i], &openouts_sk[j_indexes[i]]->alpha, values_new[i], &tmpMu, &tmpH);
        secp256k1_ge_set_gej(&out->C[i], &tmpj[i]);
        secp256k1_pedersen_commitment_save(&out->C_bin[i], &out->C[i]);

        secp256k1_scalar_get_b32(blinds[i], &openouts_sk[j_indexes[i]]->alpha);
        blind_ptr[i] = blinds[i];
        range_proof_values[i] = values_new[i];
    }

    if (alpha_c != NULL) {
        secp256k1_scalar_set_b32(&tmp, alpha_c, NULL);
        secp256k1_pedersen_ecmult(&tmpj[p_size], &tmp, withheld_c, &tmpMu, &tmpH);
        secp256k1_ge_set_gej(&out->C[p_size], &tmpj[p_size]);
        secp256k1_pedersen_commitment_save(&out->C_bin[p_size], &out->C[p_size]);

        memcpy(blinds[p_size], alpha_c, 32);
        blind_ptr[p_size] = blinds[p_size];
        range_proof_values[p_size] = withheld_c;
    }
    if (alpha_c_new != NULL) {
        int loc = p_size;
        if (alpha_c != NULL) {
            loc++;
        }
        secp256k1_scalar_set_b32(&tmp, alpha_c_new, NULL);
        secp256k1_pedersen_ecmult(&tmpj[loc], &tmp, withheld_c_new, &tmpMu, &tmpH);
        secp256k1_ge_set_gej(&out->C[loc], &tmpj[loc]);
        secp256k1_pedersen_commitment_save(&out->C_bin[loc], &out->C[loc]);

        memcpy(blinds[loc], alpha_c_new, 32);
        blind_ptr[loc] = blinds[loc];
        range_proof_values[loc] = withheld_c_new;
    }

    // range proofs
    int ret = secp256k1_bulletproof_rangeproof_prove(
            ctx,
            nctx->scratch,
            nctx->bp_gens,
            out->pi_range, &out->pi_len,
            NULL, NULL, NULL,
            range_proof_values, NULL, blind_ptr, NULL, out->range_proof_size, &nctx->mu, 64, nonce, NULL, NULL, 0, NULL);
    if (!ret) {
        free(out->C);
        free(out->C_bin);
        return 0;
    }

    out->A = (secp256k1_ge *) malloc(p_size * sizeof(secp256k1_ge));
    out->B = (secp256k1_ge *) malloc(p_size * sizeof(secp256k1_ge));
    out->f = (secp256k1_scalar *) malloc(p_size * sizeof(secp256k1_scalar));
    out->z = (secp256k1_scalar *) malloc(p_size * sizeof(secp256k1_scalar));
    out->y_a = (secp256k1_scalar *) malloc(p_size * sizeof(secp256k1_scalar));
    out->y_b = (secp256k1_scalar *) malloc(p_size * sizeof(secp256k1_scalar));

    for (i = 0; i < p_size; i++) {
        // A = C^{a}
        secp256k1_ecmult_const(&tmpj[i], &out->C[i], &openouts_sk[j_indexes[i]]->a, 256);
        secp256k1_ge_set_gej(&out->A[i], &tmpj[i]);

        // A' = h^{a'}
        secp256k1_ecmult_const(&tmpj[i], &tmpH, &openouts_sk[j_indexes[i]]->a_prime, 256);
        if (i==0) {
            secp256k1_ge_set_gej(&out->A_prime, &tmpj[i]);
        } else {
            secp256k1_gej_add_ge(&tmpj[i], &tmpj[i], &out->A_prime);
            secp256k1_ge_set_gej(&out->A_prime, &tmpj[i]);
        }

        // B = G'^{b}
        secp256k1_ecmult_const(&tmpj[i], &asset[j_indexes[i]]->G, &openouts_sk[j_indexes[i]]->b, 256);
        secp256k1_ge_set_gej(&out->B[i], &tmpj[i]);

        /* y_a */
        secp256k1_sha256_initialize(&sha);
        secp256k1_sha256_write(&sha, out->C_bin[i].data, 33);
        secp256k1_ge_save(gen_buf, &out->A[i]);
        secp256k1_sha256_write(&sha, gen_buf, 33);
        secp256k1_ge_save(gen_buf, &out->B[i]);
        secp256k1_sha256_write(&sha, gen_buf, 33);
        secp256k1_sha256_finalize(&sha, buf);
        secp256k1_scalar_set_b32(&out->y_a[i], buf, NULL);

        // fl := ρkl(v′l − vl) − ya,lal
        secp256k1_scalar_set_u64(&tmp, values[i]);
        secp256k1_scalar_negate(&tmp, &tmp);
        secp256k1_scalar_set_u64(&tmp1, values_new[i]);
        secp256k1_scalar_add(&tmp, &tmp, &tmp1);
        secp256k1_scalar_mul(&tmp1, &tmp, &acc_sk[i]->k);
        secp256k1_scalar_mul(&tmp, &out->y_a[i], &openouts_sk[j_indexes[i]]->a);
        secp256k1_scalar_negate(&tmp, &tmp);
        secp256k1_scalar_add(&out->f[i], &tmp, &tmp1);

        /* y_b */
        secp256k1_sha256_initialize(&sha);
        secp256k1_sha256_write(&sha, buf, 32);
        secp256k1_scalar_get_b32(buf, &out->f[i]);
        secp256k1_sha256_finalize(&sha, buf);
        secp256k1_scalar_set_b32(&out->y_b[i], buf, NULL);

        // zl
        secp256k1_scalar_mul(&tmp1, &tmp1, &acc_sk[i]->k);
        secp256k1_scalar_mul(&tmp, &out->y_b[i], &openouts_sk[j_indexes[i]]->b);
        secp256k1_scalar_negate(&tmp, &tmp);
        secp256k1_scalar_add(&out->z[i], &tmp1, &tmp);
    }

    // create D values
    secp256k1_ecmult_const(&tmpDj, &asset[0]->V, &openouts_sk[0]->kappa, 256);
    secp256k1_ge_set_gej(&out->D, &tmpDj);
    for (i = 1; i < n_size; i++) {
        secp256k1_ecmult_const(&tmpDj, &asset[i]->V, &openouts_sk[i]->kappa, 256);
        secp256k1_gej_add_ge(&tmpDj, &tmpDj, &out->D);
        secp256k1_ge_set_gej(&out->D, &tmpDj);
    }

    return 1;
}

void secp256k1_free_nopenena_outputs(nopenena_openout *out) {
    free(out->C);
    free(out->C_bin);
    free(out->A);
    free(out->B);
    free(out->f);
    free(out->z);
    free(out->y_a);
    free(out->y_b);
}

void secp256k1_free_nopenena_openout_proof(nopenena_openout_proof *proof) {
    free(proof->pi_zero);
}

void print_ints(int p_size, int *current_index) {
    for (int i = 0; i < p_size; i++) {
        printf("%d ", current_index[i]);
    }
    printf("\n");
}

int secp256k1_prove_nopenena_outputs2(
        const secp256k1_context* ctx,
        nopenena_context *nctx,
        nopenena_openout_proof *proof,
        nopenena_openout *out,
        nopenena_account_asset *asset[],
        uint64_t *values_new,
        uint64_t *values,
        nopenena_account_sk *acc_sk[],
        nopenena_openout_sk *openouts_sk[],
        int *j_indexes,
        int p_size,
        int n_size,
        nopenena_update_proof *uproofs[],
        unsigned char *e_bin,
        int group_model) {

    ARG_CHECK(ctx != NULL);
    ARG_CHECK(nctx != NULL);
    ARG_CHECK(proof != NULL);
    ARG_CHECK(out != NULL);
    ARG_CHECK(asset != NULL);
    ARG_CHECK(values_new != NULL);
    ARG_CHECK(acc_sk != NULL);
    ARG_CHECK(openouts_sk != NULL);
    ARG_CHECK(uproofs != NULL);
    ARG_CHECK(nctx->maximum_p_size >= p_size);
    ARG_CHECK(j_indexes[p_size - 1] <= n_size);

    secp256k1_sha256 sha;
    secp256k1_ge tmpG;
    secp256k1_ge tmpG2;
    secp256k1_ge tmpMu;
    secp256k1_ge tmpH;
    secp256k1_ge V_bar;
    secp256k1_ge P;
    secp256k1_gej tmpj;
    secp256k1_gej tmpj1;
    unsigned char buf[32];
    unsigned char gen_buf[33];
    unsigned char gen_buf0[33];
    int i;
    secp256k1_scalar y;
    secp256k1_scalar x[n_size];
    secp256k1_scalar y_a_prime;
    secp256k1_scalar tmp;
    secp256k1_scalar tmp1;
    secp256k1_scalar e;
    secp256k1_scalar xi[p_size];
    int overflow;

    secp256k1_scalar_set_b32(&e, e_bin, &overflow);
    if (overflow) {
        return 0;
    }

    /* h, mu */
    secp256k1_generator_load(&tmpH, &nctx->h);
    secp256k1_generator_load(&tmpMu, &nctx->mu);

    for (i = 0; i < n_size; i++) {
        secp256k1_scalar_set_b32(&x[i], uproofs[i]->x, NULL);
    }

    // create bar{V}
    secp256k1_ecmult_const(&tmpj, &asset[0]->V, &uproofs[0]->s3, 256);
    secp256k1_gej_add_ge(&tmpj, &tmpj, &out->D);
    secp256k1_ge_set_gej(&V_bar, &tmpj);
    for (i = 1; i < n_size; i++) {
        secp256k1_ecmult_const(&tmpj, &asset[i]->V, &uproofs[i]->s3, 256);
        secp256k1_gej_add_ge(&tmpj, &tmpj, &V_bar);
        secp256k1_ge_set_gej(&V_bar, &tmpj);
    }

    // E
    secp256k1_ecmult_const(&tmpj, &V_bar, &e, 256);
    secp256k1_ge_set_gej(&proof->E, &tmpj);

    /* y */
    secp256k1_sha256_initialize(&sha);
    secp256k1_ge_save(gen_buf, &proof->E);
    secp256k1_sha256_write(&sha, gen_buf, 33);
    for (i = 0; i < n_size; i++) {
        secp256k1_scalar_get_b32(gen_buf, &uproofs[i]->s3);
        secp256k1_sha256_write(&sha, gen_buf, 32);
        secp256k1_sha256_write(&sha, uproofs[i]->x, 32);
    }
    secp256k1_sha256_finalize(&sha, buf);
    secp256k1_scalar_set_b32(&y, buf, NULL);

    // \beta := \rho - ye
    secp256k1_scalar_mul(&tmp, &y, &e);
    secp256k1_scalar_negate(&tmp, &tmp);
    secp256k1_scalar_add(&proof->beta, &tmp, &openouts_sk[0]->rho);

    // \bar{P} := \bar{V}^{\beta}E^{y}
    secp256k1_ecmult_const(&tmpj1, &proof->E, &y, 256);
    secp256k1_ge_set_gej(&P, &tmpj1);
    secp256k1_ecmult_const(&tmpj1, &V_bar, &proof->beta, 256);
    secp256k1_gej_add_ge(&tmpj1, &tmpj1, &P);
    secp256k1_ge_set_gej(&P, &tmpj1);

    /* y_prime */
    secp256k1_sha256_initialize(&sha);
    secp256k1_ge_save(gen_buf, &P);
    secp256k1_sha256_write(&sha, gen_buf, 33);
    secp256k1_sha256_finalize(&sha, buf);
    secp256k1_scalar_set_b32(&y_a_prime, buf, NULL);

    // xi
    for (i = 0; i < p_size; i++) {
        secp256k1_scalar_set_u64(&tmp, values[i]);
        secp256k1_scalar_negate(&tmp, &tmp);
        secp256k1_scalar_set_u64(&tmp1, values_new[i]);
        secp256k1_scalar_add(&tmp, &tmp, &tmp1);
        secp256k1_scalar_mul(&tmp1, &tmp, &acc_sk[i]->k);
        secp256k1_scalar_mul(&tmp1, &tmp1, &openouts_sk[j_indexes[i]]->rho);
        secp256k1_scalar_mul(&xi[i], &tmp1, &openouts_sk[j_indexes[i]]->alpha);
        secp256k1_scalar_mul(&xi[i], &xi[i], &x[j_indexes[i]]);

        secp256k1_scalar_mul(&tmp, &y_a_prime, &openouts_sk[j_indexes[i]]->a_prime);
        secp256k1_scalar_add(&xi[i], &xi[i], &tmp);
    }

    int N;
    int j;
    int *combinations = NULL;
    if (group_model == NOPENI_PLAIN) {
        N = n_size / p_size;
        j = j_indexes[0]/p_size;
        //printf("N:%d, j:%d, log_N:%d\n", N, j, get_ceiled_log_2(N));
    } else if (group_model == NOPENI_ALL) {
        N = secp256k1_nopenena_get_total_combo(p_size, n_size);
        combinations = (int*) malloc(p_size * N * sizeof(int));
        if (!secp256k1_nopenena_get_allcombinations(combinations, p_size, n_size, N)) {
            printf("wrong total combinations\n");
            return 0;
        }
        // get correct j
        j = secp256k1_nopenena_get_j_of_all(j_indexes, combinations, p_size, N);
        if (j == -1) {
            printf("could not find j_index\n");
            return 0;
        }
        //printf("N:%d, j:%d, log_N:%d\n", N, j, get_ceiled_log_2(N));
    } else {
        printf("unknown group model\n");
        return 0;
    }

    secp256k1_pedersen_commitment H[N];
    secp256k1_ge excess;
    secp256k1_scalar xi_sum;
    secp256k1_scalar_add(&xi_sum, &xi[0], &xi[1]);
    for (i = 2; i < p_size; i++) {
        secp256k1_scalar_add(&xi_sum, &xi_sum, &xi[i]);
    }
    secp256k1_scalar_negate(&xi_sum, &xi_sum);
    secp256k1_ecmult_const(&tmpj, &tmpH, &xi_sum, 256);
    secp256k1_ge_set_gej(&excess, &tmpj);

    secp256k1_ge C_bar[p_size];
    for (i = 0; i < p_size; i++) {
        secp256k1_ecmult_const(&tmpj, &out->B[i], &out->y_b[i], 256);
        secp256k1_ge_set_gej(&C_bar[i], &tmpj);
        secp256k1_ecmult_const(&tmpj, &out->C[i], &out->f[i], 256);
        secp256k1_ge_set_gej(&tmpG2, &tmpj);
        secp256k1_ecmult_const(&tmpj, &out->A[i], &out->y_a[i], 256);
        secp256k1_gej_add_ge(&tmpj, &tmpj, &tmpG2);
        secp256k1_gej_add_ge(&tmpj, &tmpj, &C_bar[i]);
        secp256k1_ge_set_gej(&C_bar[i], &tmpj);
    }
    secp256k1_ge AY;
    secp256k1_ecmult_const(&tmpj, &out->A_prime, &y_a_prime, 256);
    secp256k1_ge_set_gej(&AY, &tmpj);

    // create H
    for (int n = 0; n < N; n++) {
        for (i = 0; i < p_size; i++) {
            int t = p_size * n + i;
            if (group_model == NOPENI_ALL) {
                t = combinations[t];
            }
            secp256k1_ecmult_const(&tmpj, &asset[t]->G, &out->z[i], 256);
            secp256k1_gej_add_ge(&tmpj, &tmpj, &C_bar[i]);
            secp256k1_ge_set_gej(&tmpG2, &tmpj);

            secp256k1_ecmult_const(&tmpj1, &tmpG2, &x[t], 256);

            if (i == 0) {
                secp256k1_gej_add_ge(&tmpj1, &tmpj1, &AY);
                secp256k1_ge_set_gej(&tmpG, &tmpj1);
            } else {
                secp256k1_gej_add_ge(&tmpj1, &tmpj1, &tmpG);
                secp256k1_ge_set_gej(&tmpG, &tmpj1);
            }
        }
        secp256k1_ge_neg(&tmpG, &tmpG);
        secp256k1_gej_set_ge(&tmpj1, &tmpG);
        secp256k1_gej_add_ge(&tmpj1, &tmpj1, &P);
        secp256k1_ge_set_gej(&tmpG, &tmpj1);
        secp256k1_pedersen_commitment_save(&H[n], &tmpG);

        if (n == j) {
            secp256k1_ge_save(gen_buf0, &excess);
            if (memcmp(H[n].data, gen_buf0, 33) != 0) {
                printf("wrong excess\n");
                return 0;
            }
        }
    }
    if (group_model == NOPENI_ALL) {
        free(combinations);
    }

    unsigned char blind[32];
    secp256k1_scalar_get_b32(blind, &xi_sum);
    int log_N = get_ceiled_log_2(N);
    proof->pi_zero_len = (32 * (log_N * 3 + 1) + 33 * log_N * 4);
    proof->pi_zero = (uint8_t *) malloc(proof->pi_zero_len * sizeof(uint8_t));
    if (!secp256k1_unlinked_logarithmic_zero_com_prove(ctx, proof->pi_zero, j, blind, N, log_N,
                                                       buf, H, &nctx->mu, &nctx->h)) {
        printf("failed zero com prove\n");
        return 0;
    }
    return 1;
}

int secp256k1_verify_nopenena_outputs2(
        const secp256k1_context* ctx,
        nopenena_context *nctx,
        nopenena_openout_proof *proof,
        nopenena_openout *out,
        nopenena_account_asset *asset[],
        int p_size,
        int n_size,
        nopenena_update_proof *uproofs[],
        int group_model) {
    ARG_CHECK(ctx != NULL);
    ARG_CHECK(nctx != NULL);
    ARG_CHECK(proof != NULL);
    ARG_CHECK(out != NULL);
    ARG_CHECK(asset != NULL);
    ARG_CHECK(uproofs != NULL);
    ARG_CHECK(nctx->maximum_p_size >= p_size);

    secp256k1_sha256 sha;
    secp256k1_ge tmpG;
    secp256k1_ge tmpG2;
    secp256k1_ge tmpMu;
    secp256k1_ge tmpH;
    secp256k1_ge V_bar;
    secp256k1_ge P;
    secp256k1_gej tmpj;
    secp256k1_gej tmpj1;
    unsigned char buf[32];
    unsigned char gen_buf[33];
    int i;
    secp256k1_scalar y;
    secp256k1_scalar x[n_size];
    secp256k1_scalar y_a_prime;


    /* h, mu */
    secp256k1_generator_load(&tmpH, &nctx->h);
    secp256k1_generator_load(&tmpMu, &nctx->mu);

    if (!secp256k1_bulletproof_rangeproof_verify(
            ctx,
            nctx->scratch,
            nctx->bp_gens,
            out->pi_range, out->pi_len,
            NULL, out->C_bin, p_size, 64, &nctx->mu, NULL, 0)) {
        printf("wrong range\n");
        return 0;
    }

    for (i = 0; i < n_size; i++) {
        secp256k1_scalar_set_b32(&x[i], uproofs[i]->x, NULL);
    }
    // create bar{V}
    secp256k1_ecmult_const(&tmpj, &asset[0]->V, &uproofs[0]->s3, 256);
    secp256k1_gej_add_ge(&tmpj, &tmpj, &out->D);
    secp256k1_ge_set_gej(&V_bar, &tmpj);
    for (i = 1; i < n_size; i++) {
        secp256k1_ecmult_const(&tmpj, &asset[i]->V, &uproofs[i]->s3, 256);
        secp256k1_gej_add_ge(&tmpj, &tmpj, &V_bar);
        secp256k1_ge_set_gej(&V_bar, &tmpj);
    }

    /* y */
    secp256k1_sha256_initialize(&sha);
    secp256k1_ge_save(gen_buf, &proof->E);
    secp256k1_sha256_write(&sha, gen_buf, 33);
    for (i = 0; i < n_size; i++) {
        secp256k1_scalar_get_b32(gen_buf, &uproofs[i]->s3);
        secp256k1_sha256_write(&sha, gen_buf, 32);
        secp256k1_sha256_write(&sha, uproofs[i]->x, 32);
    }
    secp256k1_sha256_finalize(&sha, buf);
    secp256k1_scalar_set_b32(&y, buf, NULL);

    // \bar{P} := \bar{V}^{\beta}E^{y}
    secp256k1_ecmult_const(&tmpj1, &proof->E, &y, 256);
    secp256k1_ge_set_gej(&P, &tmpj1);
    secp256k1_ecmult_const(&tmpj1, &V_bar, &proof->beta, 256);
    secp256k1_gej_add_ge(&tmpj1, &tmpj1, &P);
    secp256k1_ge_set_gej(&P, &tmpj1);

    /* y_prime */
    secp256k1_sha256_initialize(&sha);
    secp256k1_ge_save(gen_buf, &P);
    secp256k1_sha256_write(&sha, gen_buf, 33);
    secp256k1_sha256_finalize(&sha, buf);
    secp256k1_scalar_set_b32(&y_a_prime, buf, NULL);

    int N;
    int *combinations = NULL;
    if (group_model == NOPENI_PLAIN) {
        N = n_size / p_size;
        //printf("N:%d, log_N:%d\n", N, get_ceiled_log_2(N));
    } else if (group_model == NOPENI_ALL) {
        N = secp256k1_nopenena_get_total_combo(p_size, n_size);
        combinations = (int*) malloc(p_size * N * sizeof(int));
        if (!secp256k1_nopenena_get_allcombinations(combinations, p_size, n_size, N)) {
            printf("wrong total combinations\n");
            return 0;
        }
        //printf("N:%d, log_N:%d\n", N, get_ceiled_log_2(N));
    } else {
        printf("unknown group model\n");
        return 0;
    }

    secp256k1_pedersen_commitment H[N];

    secp256k1_ge C_bar[p_size];
    for (i = 0; i < p_size; i++) {
        secp256k1_ecmult_const(&tmpj, &out->B[i], &out->y_b[i], 256);
        secp256k1_ge_set_gej(&C_bar[i], &tmpj);
        secp256k1_ecmult_const(&tmpj, &out->C[i], &out->f[i], 256);
        secp256k1_ge_set_gej(&tmpG2, &tmpj);
        secp256k1_ecmult_const(&tmpj, &out->A[i], &out->y_a[i], 256);
        secp256k1_gej_add_ge(&tmpj, &tmpj, &tmpG2);
        secp256k1_gej_add_ge(&tmpj, &tmpj, &C_bar[i]);
        secp256k1_ge_set_gej(&C_bar[i], &tmpj);
    }
    secp256k1_ge AY;
    secp256k1_ecmult_const(&tmpj, &out->A_prime, &y_a_prime, 256);
    secp256k1_ge_set_gej(&AY, &tmpj);

    // create H
    for (int n = 0; n < N; n++) {
        for (i = 0; i < p_size; i++) {
            int t = p_size * n + i;
            if (group_model == NOPENI_ALL) {
                t = combinations[t];
            }
            secp256k1_ecmult_const(&tmpj, &asset[t]->G, &out->z[i], 256);
            secp256k1_gej_add_ge(&tmpj, &tmpj, &C_bar[i]);
            secp256k1_ge_set_gej(&tmpG2, &tmpj);

            secp256k1_ecmult_const(&tmpj1, &tmpG2, &x[t], 256);

            if (i == 0) {
                secp256k1_gej_add_ge(&tmpj1, &tmpj1, &AY);
                secp256k1_ge_set_gej(&tmpG, &tmpj1);
            } else {
                secp256k1_gej_add_ge(&tmpj1, &tmpj1, &tmpG);
                secp256k1_ge_set_gej(&tmpG, &tmpj1);
            }
        }
        secp256k1_ge_neg(&tmpG, &tmpG);
        secp256k1_gej_set_ge(&tmpj1, &tmpG);
        secp256k1_gej_add_ge(&tmpj1, &tmpj1, &P);
        secp256k1_ge_set_gej(&tmpG, &tmpj1);
        secp256k1_pedersen_commitment_save(&H[n], &tmpG);
    }
    if (group_model == NOPENI_ALL) {
        free(combinations);
    }

    int log_N = get_ceiled_log_2(N);
    return secp256k1_unlinked_logarithmic_zero_com_verify(ctx, proof->pi_zero, N, log_N,
                                                          buf, H, &nctx->mu, &nctx->h);
}

int secp256k1_prove_nopenena_outputs(
        const secp256k1_context* ctx,
        nopenena_context *nctx,
        nopenena_openout_proof *proof,
        nopenena_openout *out,
        nopenena_account_asset *asset[],
        uint64_t *values_new,
        uint64_t *values,
        nopenena_account_sk *acc_sk[],
        nopenena_openout_sk *openouts_sk[],
        int *j_indexes,
        int p_size,
        int n_size,
        nopenena_update_proof *uproofs[],
        unsigned char *e_bin,
        int group_model,
        int zero_com_flag) {

    ARG_CHECK(ctx != NULL);
    ARG_CHECK(nctx != NULL);
    ARG_CHECK(proof != NULL);
    ARG_CHECK(out != NULL);
    ARG_CHECK(asset != NULL);
    ARG_CHECK(values_new != NULL);
    ARG_CHECK(acc_sk != NULL);
    ARG_CHECK(openouts_sk != NULL);
    ARG_CHECK(uproofs != NULL);
    ARG_CHECK(nctx->maximum_p_size >= p_size);
    ARG_CHECK(j_indexes[p_size - 1] <= n_size);

    secp256k1_sha256 sha;
    secp256k1_ge tmpG;
    secp256k1_ge tmpG2;
    secp256k1_ge tmpMu;
    secp256k1_ge tmpH;
    secp256k1_ge V_bar;
    secp256k1_gej tmpj;
    secp256k1_gej tmpj1;
    unsigned char buf[32];
    unsigned char gen_buf[33];
    unsigned char gen_buf0[33];
    int i;
    secp256k1_scalar y;
    secp256k1_scalar x[n_size];
    secp256k1_scalar tmp;
    secp256k1_scalar tmp1;
    secp256k1_scalar e;
    secp256k1_scalar xi[p_size];
    int overflow;

    secp256k1_scalar_set_b32(&e, e_bin, &overflow);
    if (overflow) {
        return 0;
    }

    /* h, mu */
    secp256k1_generator_load(&tmpH, &nctx->h);
    secp256k1_generator_load(&tmpMu, &nctx->mu);

    for (i = 0; i < n_size; i++) {
        secp256k1_scalar_set_b32(&x[i], uproofs[i]->x, NULL);
    }

    // create bar{V}
    secp256k1_ecmult_const(&tmpj, &asset[0]->V, &uproofs[0]->s3, 256);
    secp256k1_gej_add_ge(&tmpj, &tmpj, &out->D);
    secp256k1_ge_set_gej(&V_bar, &tmpj);
    for (i = 1; i < n_size; i++) {
        secp256k1_ecmult_const(&tmpj, &asset[i]->V, &uproofs[i]->s3, 256);
        secp256k1_gej_add_ge(&tmpj, &tmpj, &V_bar);
        secp256k1_ge_set_gej(&V_bar, &tmpj);
    }

    /* y */
    secp256k1_sha256_initialize(&sha);
    secp256k1_ge_save(gen_buf, &V_bar);
    secp256k1_sha256_write(&sha, gen_buf, 33);
    for (i = 0; i < n_size; i++) {
        secp256k1_scalar_get_b32(gen_buf, &uproofs[i]->s3);
        secp256k1_sha256_write(&sha, gen_buf, 32);
        secp256k1_sha256_write(&sha, uproofs[i]->x, 32);
    }
    secp256k1_sha256_finalize(&sha, buf);
    secp256k1_scalar_set_b32(&y, buf, NULL);

    // xi
    for (i = 0; i < p_size; i++) {
        secp256k1_scalar_set_u64(&tmp, values[i]);
        secp256k1_scalar_negate(&tmp, &tmp);
        secp256k1_scalar_set_u64(&tmp1, values_new[i]);
        secp256k1_scalar_add(&tmp, &tmp, &tmp1);
        secp256k1_scalar_mul(&tmp1, &tmp, &acc_sk[i]->k);
        secp256k1_scalar_mul(&xi[i], &tmp1, &openouts_sk[j_indexes[i]]->alpha);
        secp256k1_scalar_mul(&xi[i], &xi[i], &x[j_indexes[i]]);

        secp256k1_scalar_mul(&tmp, &y, &openouts_sk[j_indexes[i]]->a_prime);
        secp256k1_scalar_add(&xi[i], &xi[i], &tmp);
    }

    int N;
    int j;
    int *combinations = NULL;
    if (group_model == NOPENI_PLAIN) {
        N = n_size / p_size;
        j = j_indexes[0]/p_size;
        //printf("N:%d, j:%d, log_N:%d\n", N, j, get_ceiled_log_2(N));
    } else if (group_model == NOPENI_ALL) {
        N = secp256k1_nopenena_get_total_combo(p_size, n_size);
        combinations = (int*) malloc(p_size * N * sizeof(int));
        if (!secp256k1_nopenena_get_allcombinations(combinations, p_size, n_size, N)) {
            printf("wrong total combinations\n");
            return 0;
        }
        // get correct j
        j = secp256k1_nopenena_get_j_of_all(j_indexes, combinations, p_size, N);
        if (j == -1) {
            printf("could not find j_index\n");
            return 0;
        }
        //printf("N:%d, j:%d, log_N:%d\n", N, j, get_ceiled_log_2(N));
    } else {
        printf("unknown group model\n");
        return 0;
    }

    secp256k1_pedersen_commitment H[N];
    secp256k1_ge excess;
    secp256k1_scalar xi_sum;
    //secp256k1_scalar_add(&xi_sum, &xi[0], &xi[1]);
    uint8_t x1buf[32];
    secp256k1_scalar_get_b32(x1buf, &xi[0]);
    secp256k1_scalar_set_b32(&xi_sum, x1buf, NULL);
    for (i = 1; i < p_size; i++) {
        secp256k1_scalar_add(&xi_sum, &xi_sum, &xi[i]);
    }
    secp256k1_scalar_negate(&xi_sum, &xi_sum);
    secp256k1_ecmult_const(&tmpj, &tmpH, &xi_sum, 256);
    secp256k1_ge_set_gej(&excess, &tmpj);

    secp256k1_ge C_bar[p_size];
    for (i = 0; i < p_size; i++) {
        secp256k1_ecmult_const(&tmpj, &out->B[i], &out->y_b[i], 256);
        secp256k1_ge_set_gej(&C_bar[i], &tmpj);
        secp256k1_ecmult_const(&tmpj, &out->C[i], &out->f[i], 256);
        secp256k1_ge_set_gej(&tmpG2, &tmpj);
        secp256k1_ecmult_const(&tmpj, &out->A[i], &out->y_a[i], 256);
        secp256k1_gej_add_ge(&tmpj, &tmpj, &tmpG2);
        secp256k1_gej_add_ge(&tmpj, &tmpj, &C_bar[i]);
        secp256k1_ge_set_gej(&C_bar[i], &tmpj);
    }
    secp256k1_ge AY;
    secp256k1_ecmult_const(&tmpj, &out->A_prime, &y, 256);
    secp256k1_ge_set_gej(&AY, &tmpj);

    // create H
    for (int n = 0; n < N; n++) {
        for (i = 0; i < p_size; i++) {
            int t = p_size * n + i;
            if (group_model == NOPENI_ALL) {
                t = combinations[t];
            }
            secp256k1_ecmult_const(&tmpj, &asset[t]->G, &out->z[i], 256);
            secp256k1_gej_add_ge(&tmpj, &tmpj, &C_bar[i]);
            secp256k1_ge_set_gej(&tmpG2, &tmpj);

            secp256k1_ecmult_const(&tmpj1, &tmpG2, &x[t], 256);

            if (i == 0) {
                secp256k1_gej_add_ge(&tmpj1, &tmpj1, &AY);
                secp256k1_ge_set_gej(&tmpG, &tmpj1);
            } else {
                secp256k1_gej_add_ge(&tmpj1, &tmpj1, &tmpG);
                secp256k1_ge_set_gej(&tmpG, &tmpj1);
            }
        }
        secp256k1_ge_neg(&tmpG, &tmpG);
        secp256k1_gej_set_ge(&tmpj1, &tmpG);
        secp256k1_gej_add_ge(&tmpj1, &tmpj1, &V_bar);
        secp256k1_ge_set_gej(&tmpG, &tmpj1);
        secp256k1_pedersen_commitment_save(&H[n], &tmpG);

        if (n == j) {
            secp256k1_ge_save(gen_buf0, &excess);
            if (memcmp(H[n].data, gen_buf0, 33) != 0) {
                printf("wrong excess %d %d %d\n", n, j, N);
                return 0;
            }
        }
    }
    if (group_model == NOPENI_ALL) {
        free(combinations);
    }

    unsigned char blind[32];
    secp256k1_scalar_get_b32(blind, &xi_sum);
    int log_N = get_ceiled_log_2(N);
    if (zero_com_flag == 0) {
        proof->pi_zero_len = (32 * (log_N * 3 + 1) + 33 * log_N * 4);
        proof->pi_zero = (uint8_t *) malloc(proof->pi_zero_len * sizeof(uint8_t));
        if (!secp256k1_unlinked_logarithmic_zero_com_prove(ctx, proof->pi_zero, j, blind, N, log_N,
                                                           buf, H, &nctx->mu, &nctx->h)) {
            printf("failed zero com prove\n");
            return 0;
        }
    } else {
        proof->pi_zero_len = secp256k1_zero_mcom_get_size(&nctx->rctx);
        proof->pi_zero = (uint8_t *) malloc(proof->pi_zero_len * sizeof(uint8_t));
        if (!secp256k1_create_zero_mcom_short(ctx,
                                              &nctx->rctx,
                                              proof->pi_zero,
                                              j,
                                              blind,
                                              buf,
                                              H, N)) {
            printf("failed zero com prove\n");
            return 0;
        }
        if (!secp256k1_verify_zero_mcom_short(ctx,
                                              &nctx->rctx,
                                              proof->pi_zero,
                                              buf,
                                              H, N)) {
            printf("failed zero com verify\n");
            return 0;
        }
    }
    return 1;
}


int secp256k1_verify_nopenena_outputs(
        const secp256k1_context* ctx,
        nopenena_context *nctx,
        nopenena_openout_proof *proof,
        nopenena_openout *out,
        nopenena_account_asset *asset[],
        int p_size,
        int n_size,
        nopenena_update_proof *uproofs[],
        int group_model,
        int zero_com_flag) {
    ARG_CHECK(ctx != NULL);
    ARG_CHECK(nctx != NULL);
    ARG_CHECK(proof != NULL);
    ARG_CHECK(out != NULL);
    ARG_CHECK(asset != NULL);
    ARG_CHECK(uproofs != NULL);
    ARG_CHECK(nctx->maximum_p_size >= p_size);

    secp256k1_sha256 sha;
    secp256k1_ge tmpG;
    secp256k1_ge tmpG2;
    secp256k1_ge tmpMu;
    secp256k1_ge tmpH;
    secp256k1_ge V_bar;
    secp256k1_gej tmpj;
    secp256k1_gej tmpj1;
    unsigned char buf[32];
    unsigned char gen_buf[33];
    int i;
    secp256k1_scalar y;
    secp256k1_scalar x[n_size];


    /* h, mu */
    secp256k1_generator_load(&tmpH, &nctx->h);
    secp256k1_generator_load(&tmpMu, &nctx->mu);


    if (!secp256k1_bulletproof_rangeproof_verify(
            ctx,
            nctx->scratch,
            nctx->bp_gens,
            out->pi_range, out->pi_len,
            NULL, out->C_bin, out->range_proof_size, 64, &nctx->mu, NULL, 0)) {
        printf("wrong range due to 2^x neq %d\n", out->range_proof_size);
        return 0;
    }

    for (i = 0; i < n_size; i++) {
        secp256k1_scalar_set_b32(&x[i], uproofs[i]->x, NULL);
    }
    // create bar{V}
    secp256k1_ecmult_const(&tmpj, &asset[0]->V, &uproofs[0]->s3, 256);
    secp256k1_gej_add_ge(&tmpj, &tmpj, &out->D);
    secp256k1_ge_set_gej(&V_bar, &tmpj);
    for (i = 1; i < n_size; i++) {
        secp256k1_ecmult_const(&tmpj, &asset[i]->V, &uproofs[i]->s3, 256);
        secp256k1_gej_add_ge(&tmpj, &tmpj, &V_bar);
        secp256k1_ge_set_gej(&V_bar, &tmpj);
    }

    /* y */
    secp256k1_sha256_initialize(&sha);
    secp256k1_ge_save(gen_buf, &V_bar);
    secp256k1_sha256_write(&sha, gen_buf, 33);
    for (i = 0; i < n_size; i++) {
        secp256k1_scalar_get_b32(gen_buf, &uproofs[i]->s3);
        secp256k1_sha256_write(&sha, gen_buf, 32);
        secp256k1_sha256_write(&sha, uproofs[i]->x, 32);
    }
    secp256k1_sha256_finalize(&sha, buf);
    secp256k1_scalar_set_b32(&y, buf, NULL);


    int N;
    int *combinations = NULL;
    if (group_model == NOPENI_PLAIN) {
        N = n_size / p_size;
        //printf("N:%d, log_N:%d\n", N, get_ceiled_log_2(N));
    } else if (group_model == NOPENI_ALL) {
        N = secp256k1_nopenena_get_total_combo(p_size, n_size);
        combinations = (int*) malloc(p_size * N * sizeof(int));
        if (!secp256k1_nopenena_get_allcombinations(combinations, p_size, n_size, N)) {
            printf("wrong total combinations\n");
            return 0;
        }
        //printf("N:%d, log_N:%d\n", N, get_ceiled_log_2(N));
    } else {
        printf("unknown group model\n");
        return 0;
    }

    secp256k1_pedersen_commitment H[N];

    secp256k1_ge C_bar[p_size];
    for (i = 0; i < p_size; i++) {
        secp256k1_ecmult_const(&tmpj, &out->B[i], &out->y_b[i], 256);
        secp256k1_ge_set_gej(&C_bar[i], &tmpj);
        secp256k1_ecmult_const(&tmpj, &out->C[i], &out->f[i], 256);
        secp256k1_ge_set_gej(&tmpG2, &tmpj);
        secp256k1_ecmult_const(&tmpj, &out->A[i], &out->y_a[i], 256);
        secp256k1_gej_add_ge(&tmpj, &tmpj, &tmpG2);
        secp256k1_gej_add_ge(&tmpj, &tmpj, &C_bar[i]);
        secp256k1_ge_set_gej(&C_bar[i], &tmpj);
    }
    secp256k1_ge AY;
    secp256k1_ecmult_const(&tmpj, &out->A_prime, &y, 256);
    secp256k1_ge_set_gej(&AY, &tmpj);

    // create H
    for (int n = 0; n < N; n++) {
        for (i = 0; i < p_size; i++) {
            int t = p_size * n + i;
            if (group_model == NOPENI_ALL) {
                t = combinations[t];
            }
            secp256k1_ecmult_const(&tmpj, &asset[t]->G, &out->z[i], 256);
            secp256k1_gej_add_ge(&tmpj, &tmpj, &C_bar[i]);
            secp256k1_ge_set_gej(&tmpG2, &tmpj);

            secp256k1_ecmult_const(&tmpj1, &tmpG2, &x[t], 256);

            if (i == 0) {
                secp256k1_gej_add_ge(&tmpj1, &tmpj1, &AY);
                secp256k1_ge_set_gej(&tmpG, &tmpj1);
            } else {
                secp256k1_gej_add_ge(&tmpj1, &tmpj1, &tmpG);
                secp256k1_ge_set_gej(&tmpG, &tmpj1);
            }
        }
        secp256k1_ge_neg(&tmpG, &tmpG);
        secp256k1_gej_set_ge(&tmpj1, &tmpG);
        secp256k1_gej_add_ge(&tmpj1, &tmpj1, &V_bar);
        secp256k1_ge_set_gej(&tmpG, &tmpj1);
        secp256k1_pedersen_commitment_save(&H[n], &tmpG);
    }
    if (group_model == NOPENI_ALL) {
        free(combinations);
    }

    int log_N = get_ceiled_log_2(N);
    if (zero_com_flag == 0) {
        return secp256k1_unlinked_logarithmic_zero_com_verify(ctx,
                                                              proof->pi_zero,
                                                              N, log_N,
                                                              buf, H,
                                                              &nctx->mu,
                                                              &nctx->h);
    }
    return secp256k1_verify_zero_mcom_short(ctx,
                                             &nctx->rctx,
                                             proof->pi_zero,
                                             buf,
                                             H, N);
}


int secp256k1_prove_nopenena_balance(
        const secp256k1_context* ctx,
        const nopenena_context *nctx,
        nopenena_balance_proof *proof,
        uint64_t f,
        uint64_t f_new,
        nopenena_account *acc,
        nopenena_account_asset *asset0[],
        nopenena_account_asset *asset1[],
        const uint64_t *values,
        const uint64_t *values_new,
        int n_size,
        unsigned char *r_bin,
        unsigned char *u_bin,
        unsigned char *u_prime_bin,
        secp256k1_pedersen_commitment *C,
        unsigned char *alpha_c,
        uint64_t withheld_c,
        secp256k1_pedersen_commitment *C_new,
        unsigned char *alpha_c_new,
        uint64_t withheld_c_new) {
    ARG_CHECK(ctx != NULL);
    ARG_CHECK(proof != NULL);
    ARG_CHECK(acc != NULL);
    ARG_CHECK(asset0 != NULL);
    ARG_CHECK(asset1 != NULL);

    int i, overflow;
    secp256k1_ge tmpG;
    secp256k1_ge tmpMu;
    secp256k1_ge tmpH;
    secp256k1_ge E;
    secp256k1_ge K;
    secp256k1_gej tmpj;
    secp256k1_scalar tmp;
    secp256k1_scalar tmp1;
    secp256k1_scalar u;
    secp256k1_scalar u_prime;
    secp256k1_scalar r;
    secp256k1_scalar y;
    unsigned char gen_buf[33];
    secp256k1_sha256 sha;

    uint64_t inputs = 0;
    uint64_t outputs = 0;
    for (i =0; i < n_size; i++) {
        inputs += values[i];
        outputs += values_new[i];
    }
    inputs += f; // reward
    outputs += f_new; // fee
    inputs += withheld_c;
    outputs += withheld_c_new;

    if (inputs != outputs) {
        printf("invalid total coins %ld %ld\n", inputs, outputs);
        return 0;
    }

    secp256k1_scalar_set_b32(&u, u_bin, &overflow);
    if (overflow) {
        return 0;
    }
    secp256k1_scalar_set_b32(&u_prime, u_prime_bin, &overflow);
    if (overflow) {
        return 0;
    }
    secp256k1_scalar_set_b32(&r, r_bin, &overflow);
    if (overflow) {
        return 0;
    }

    secp256k1_generator_load(&tmpMu, &nctx->mu);
    secp256k1_generator_load(&tmpH, &nctx->h);
    // E
    for (i = 0; i < n_size; i++) {
        secp256k1_ge_neg(&tmpG, &asset0[i]->V);
        secp256k1_gej_set_ge(&tmpj, &tmpG);
        secp256k1_gej_add_ge(&tmpj, &tmpj, &asset1[i]->V);

        if (i != 0) {
            secp256k1_gej_add_ge(&tmpj, &tmpj, &E);
            secp256k1_ge_set_gej(&E, &tmpj);
        } else {
            secp256k1_ge_set_gej(&E, &tmpj);
        }
    }
    if (C_new != NULL) {
        secp256k1_ge_load(C_new->data, &tmpG);
        //secp256k1_ge_neg(&tmpG, &tmpG);
        secp256k1_gej_set_ge(&tmpj, &tmpG);
        secp256k1_gej_add_ge(&tmpj, &tmpj, &E);
        secp256k1_ge_set_gej(&E, &tmpj);
    }
    if (C != NULL) {
        secp256k1_ge_load(C->data, &tmpG);
        secp256k1_ge_neg(&tmpG, &tmpG);
        secp256k1_gej_set_ge(&tmpj, &tmpG);
        secp256k1_gej_add_ge(&tmpj, &tmpj, &E);
        secp256k1_ge_set_gej(&E, &tmpj);
    }

    secp256k1_scalar_set_u64(&tmp, f);
    secp256k1_scalar_set_u64(&tmp1, f_new);
    secp256k1_scalar_negate(&tmp, &tmp);
    secp256k1_scalar_add(&tmp, &tmp, &tmp1);

    secp256k1_ecmult_const(&tmpj, &tmpMu, &tmp, 256);
    secp256k1_gej_add_ge(&tmpj, &tmpj, &E);
    secp256k1_ge_set_gej(&E, &tmpj);

    // K
    secp256k1_gej_set_ge(&tmpj, &acc[0].K);
    for (i = 1; i < n_size; i++) {
        secp256k1_gej_add_ge(&tmpj, &tmpj, &acc[i].K);
    }
    secp256k1_ge_set_gej(&K, &tmpj);

    // U
    secp256k1_ecmult_const(&tmpj, &K, &u, 256);
    secp256k1_ge_set_gej(&proof->U, &tmpj);
    secp256k1_ecmult_const(&tmpj, &tmpH, &u_prime, 256);
    secp256k1_gej_add_ge(&tmpj, &tmpj, &proof->U);
    secp256k1_ge_set_gej(&proof->U, &tmpj);

    /* y */
    secp256k1_sha256_initialize(&sha);
    secp256k1_ge_save(gen_buf, &E);
    secp256k1_sha256_write(&sha, gen_buf, 33);
    secp256k1_ge_save(gen_buf, &K);
    secp256k1_sha256_write(&sha, gen_buf, 33);
    secp256k1_ge_save(gen_buf, &proof->U);
    secp256k1_sha256_write(&sha, gen_buf, 33);
    secp256k1_sha256_finalize(&sha, gen_buf);
    secp256k1_scalar_set_b32(&y, gen_buf, NULL);

    secp256k1_scalar_mul(&tmp, &y, &u);
    secp256k1_scalar_negate(&tmp, &tmp);
    secp256k1_scalar_add(&proof->s_bar, &tmp, &r);

    secp256k1_scalar_mul(&tmp, &y, &u_prime);
    secp256k1_scalar_negate(&proof->s_bar_prime, &tmp);
    if (alpha_c_new != NULL) {
        secp256k1_scalar_set_b32(&tmp, alpha_c_new, NULL);
        secp256k1_scalar_add(&proof->s_bar_prime, &proof->s_bar_prime, &tmp);
    }
    if (alpha_c != NULL) {
        secp256k1_scalar_set_b32(&tmp, alpha_c, NULL);
        secp256k1_scalar_negate(&tmp, &tmp);
        secp256k1_scalar_add(&proof->s_bar_prime, &proof->s_bar_prime, &tmp);
    }

    /*secp256k1_ecmult_const(&tmpj, &proof->U, &y, 256);
    secp256k1_ge_set_gej(&tmpG, &tmpj);
    secp256k1_ecmult_const(&tmpj, &K, &proof->s_bar, 256);
    secp256k1_gej_add_ge(&tmpj, &tmpj, &tmpG);
    secp256k1_ge_set_gej(&tmpG, &tmpj);
    secp256k1_ecmult_const(&tmpj, &tmpH, &proof->s_bar_prime, 256);
    secp256k1_gej_add_ge(&tmpj, &tmpj, &tmpG);
    secp256k1_ge_set_gej(&tmpG, &tmpj);


    unsigned char gen_buf0[33];
    secp256k1_ge_save(gen_buf, &E);
    secp256k1_ge_save(gen_buf0, &tmpG);
    printf("val %d\n", memcmp(gen_buf, gen_buf0, 33) == 0);*/

    return 1;
}

int secp256k1_verify_nopenena_balance(
        const secp256k1_context* ctx,
        const nopenena_context *nctx,
        nopenena_balance_proof *proof,
        uint64_t f,
        uint64_t f_new,
        nopenena_account *acc,
        nopenena_account_asset *asset0[],
        nopenena_account_asset *asset1[],
        int n_size,
        secp256k1_pedersen_commitment *C,
        secp256k1_pedersen_commitment *C_new) {
    ARG_CHECK(ctx != NULL);
    ARG_CHECK(proof != NULL);
    ARG_CHECK(acc != NULL);
    ARG_CHECK(asset0 != NULL);
    ARG_CHECK(asset1 != NULL);


    int i;
    secp256k1_ge tmpG;
    secp256k1_ge tmpMu;
    secp256k1_ge tmpH;
    secp256k1_ge E;
    secp256k1_ge K;
    secp256k1_gej tmpj;
    secp256k1_scalar y;
    secp256k1_scalar tmp;
    secp256k1_scalar tmp1;
    unsigned char gen_buf[33];
    unsigned char gen_buf0[33];
    secp256k1_sha256 sha;


    secp256k1_generator_load(&tmpMu, &nctx->mu);
    secp256k1_generator_load(&tmpH, &nctx->h);
    secp256k1_generator_load(&tmpMu, &nctx->mu);
    secp256k1_generator_load(&tmpH, &nctx->h);
    // E
    for (i = 0; i < n_size; i++) {
        secp256k1_ge_neg(&tmpG, &asset0[i]->V);
        secp256k1_gej_set_ge(&tmpj, &tmpG);
        secp256k1_gej_add_ge(&tmpj, &tmpj, &asset1[i]->V);

        if (i != 0) {
            secp256k1_gej_add_ge(&tmpj, &tmpj, &E);
            secp256k1_ge_set_gej(&E, &tmpj);
        } else {
            secp256k1_ge_set_gej(&E, &tmpj);
        }
    }
    if (C_new != NULL) {
        secp256k1_ge_load(C_new->data, &tmpG);
        //secp256k1_ge_neg(&tmpG, &tmpG);
        secp256k1_gej_set_ge(&tmpj, &tmpG);
        secp256k1_gej_add_ge(&tmpj, &tmpj, &E);
        secp256k1_ge_set_gej(&E, &tmpj);
    }
    if (C != NULL) {
        secp256k1_ge_load(C->data, &tmpG);
        secp256k1_ge_neg(&tmpG, &tmpG);
        secp256k1_gej_set_ge(&tmpj, &tmpG);
        secp256k1_gej_add_ge(&tmpj, &tmpj, &E);
        secp256k1_ge_set_gej(&E, &tmpj);
    }

    secp256k1_scalar_set_u64(&tmp, f);
    secp256k1_scalar_set_u64(&tmp1, f_new);
    secp256k1_scalar_negate(&tmp, &tmp);
    secp256k1_scalar_add(&tmp, &tmp, &tmp1);

    secp256k1_ecmult_const(&tmpj, &tmpMu, &tmp, 256);
    secp256k1_gej_add_ge(&tmpj, &tmpj, &E);
    secp256k1_ge_set_gej(&E, &tmpj);

    // K
    secp256k1_gej_set_ge(&tmpj, &acc[0].K);
    for (i = 1; i < n_size; i++) {
        secp256k1_gej_add_ge(&tmpj, &tmpj, &acc[i].K);
    }
    secp256k1_ge_set_gej(&K, &tmpj);

    /* y */
    secp256k1_sha256_initialize(&sha);
    secp256k1_ge_save(gen_buf, &E);
    secp256k1_sha256_write(&sha, gen_buf, 33);
    secp256k1_ge_save(gen_buf, &K);
    secp256k1_sha256_write(&sha, gen_buf, 33);
    secp256k1_ge_save(gen_buf, &proof->U);
    secp256k1_sha256_write(&sha, gen_buf, 33);
    secp256k1_sha256_finalize(&sha, gen_buf);
    secp256k1_scalar_set_b32(&y, gen_buf, NULL);

    secp256k1_ecmult_const(&tmpj, &proof->U, &y, 256);
    secp256k1_ge_set_gej(&tmpG, &tmpj);
    secp256k1_ecmult_const(&tmpj, &K, &proof->s_bar, 256);
    secp256k1_gej_add_ge(&tmpj, &tmpj, &tmpG);
    secp256k1_ge_set_gej(&tmpG, &tmpj);
    secp256k1_ecmult_const(&tmpj, &tmpH, &proof->s_bar_prime, 256);
    secp256k1_gej_add_ge(&tmpj, &tmpj, &tmpG);
    secp256k1_ge_set_gej(&tmpG, &tmpj);

    secp256k1_ge_save(gen_buf, &E);
    secp256k1_ge_save(gen_buf0, &tmpG);

    return memcmp(gen_buf, gen_buf0, 33) == 0;
}

#endif
