#ifndef _SECP256K1_NOPENI_
# define _SECP256K1_NOPENI_

# include "secp256k1.h"
#include "secp256k1_ringcip.h"

# ifdef __cplusplus
extern "C" {
# endif

#define NOPENI_PLAIN    0 // p_size * M = acc_size and there are M groups
#define NOPENI_ALL    1 // all possible groups

typedef struct nopenena_context_struct{
    secp256k1_generator g;
    secp256k1_generator mu;
    secp256k1_generator h;
    ringcip_context rctx;
    secp256k1_bulletproof_generators *bp_gens;
    secp256k1_context *none;
    secp256k1_scratch *scratch;
    int maximum_p_size;
} nopenena_context;

typedef struct nopenena_account_struct{
    secp256k1_ge R;
    secp256k1_ge K;
} nopenena_account;

typedef struct nopenena_account_asset_struct{
    secp256k1_ge G;
    secp256k1_ge V;
} nopenena_account_asset;

typedef struct nopenena_account_sk_struct{
    secp256k1_scalar r;
    secp256k1_scalar k;
    secp256k1_scalar gamma;
} nopenena_account_sk;

typedef struct nopenena_update_proof_struct{
    unsigned char x[32];
    secp256k1_scalar s1;
    secp256k1_scalar s2;
    secp256k1_scalar s3;
} nopenena_update_proof;

typedef struct nopenena_openout_struct{
    secp256k1_ge *C;
    secp256k1_pedersen_commitment *C_bin;
    secp256k1_ge *A;
    secp256k1_ge *B;
    secp256k1_scalar *f;
    secp256k1_scalar *z;
    secp256k1_ge D;
    secp256k1_ge A_prime;
    unsigned char pi_range[2000];
    size_t pi_len;
    unsigned char pi_contract[2000];
    size_t pi_contract_len;
    secp256k1_scalar *y_a; // will not be included in the proof
    secp256k1_scalar *y_b; // will not be included in the proof
    int p_size;
    int range_proof_size;
} nopenena_openout;

typedef struct nopenena_openout_sk_struct{
    secp256k1_scalar alpha;
    secp256k1_scalar a;
    secp256k1_scalar a_prime;
    secp256k1_scalar b;
    secp256k1_scalar kappa;
    secp256k1_scalar rho; // shared between the participants
} nopenena_openout_sk;

typedef struct nopenena_openout_proof_struct{
    secp256k1_ge E;
    secp256k1_scalar beta;
    unsigned char *pi_zero;
    int pi_zero_len;
} nopenena_openout_proof;

typedef struct nopenena_balance_proof_struct{
    secp256k1_ge U;
    secp256k1_scalar s_bar;
    secp256k1_scalar s_bar_prime;
} nopenena_balance_proof;


typedef struct {
    unsigned char data[32];
} secp256k1_hash;

typedef struct {
    unsigned char data[64];
} secp256k1_identity_pf;

typedef struct {
    unsigned char *data;
    unsigned len;
} secp256k1_unlinked_identity_pf;


typedef struct nopenena_escrow_zk_struct{
    secp256k1_pedersen_commitment E;
    secp256k1_pedersen_commitment CS;
    secp256k1_pedersen_commitment CB;
} nopenena_escrow_zk;



/**
 * Create an identity from the digest
 * @param ctx  - context object
 * @param id - output
 * @param blind
 * @param digest
 * @param value_gen
 * @param blind_gen
 * @return 1 - successful, 0 - unsuccessful (most probably because of invalid key, e.g., a blind is zero
 */
SECP256K1_API SECP256K1_WARN_UNUSED_RESULT int secp256k1_commitment_create(
        const secp256k1_context* ctx,
        secp256k1_pedersen_commitment *com,
        const unsigned char *blind,
        secp256k1_hash digest,
        const secp256k1_generator *value_gen,
        const secp256k1_generator *blind_gen
) SECP256K1_ARG_NONNULL(1) SECP256K1_ARG_NONNULL(2) SECP256K1_ARG_NONNULL(3) SECP256K1_ARG_NONNULL(5) SECP256K1_ARG_NONNULL(6);


/**
 * Create a BDPoE for an identity
 * @param ctx  - context object
 * @param proof - out
 * @param index - corresponding identity index
 * @param blind - secret key
 * @param nonce
 * @param digest
 * @param challenge
 * @param N - number of commitments
 * @param n - rounded up log(N) (can be more than log(N))
 * @param ids - all identities
 * @param value_gen
 * @param blind_gen
 * @return 1 - successful, 0 - unsuccessful (most probably because of invalid key, e.g., a blind is zero
 */
SECP256K1_API SECP256K1_WARN_UNUSED_RESULT int secp256k1_unlinked_logarithmic_identity_prove(
        const secp256k1_context* ctx,
        secp256k1_unlinked_identity_pf *proof,
        int index,
        const unsigned char *blind,
        const unsigned char *digest,
        const unsigned char *challenge,
        int N,
        int n,
        secp256k1_pedersen_commitment *coms,
        const secp256k1_generator *value_gen,
        const secp256k1_generator *blind_gen
) SECP256K1_ARG_NONNULL(1) SECP256K1_ARG_NONNULL(2) SECP256K1_ARG_NONNULL(4) SECP256K1_ARG_NONNULL(5) SECP256K1_ARG_NONNULL(6)
SECP256K1_ARG_NONNULL(10) SECP256K1_ARG_NONNULL(11) SECP256K1_ARG_NONNULL(12);

/**
 * Verify a DBPoE of an identity
 * @param ctx  - context object
 * @param proof - DBPoE
 * @param digest
 * @param challenge
 * @param others
 * @param N - number of commitments
 * @param n - rounded up log(N) (can be more than log(N))
 * @param value_gen
 * @param blind_gen
 * @return 1 - valid, 0 - invalid
 */
SECP256K1_API SECP256K1_WARN_UNUSED_RESULT int secp256k1_unlinked_logarithmic_identity_verify(
        const secp256k1_context* ctx,
        secp256k1_unlinked_identity_pf *proof,
        const unsigned char *digest,
        const unsigned char *challenge,
        secp256k1_pedersen_commitment *others,
        int N,
        int n,
        const secp256k1_generator *value_gen,
        const secp256k1_generator *blind_gen
) SECP256K1_ARG_NONNULL(1) SECP256K1_ARG_NONNULL(2) SECP256K1_ARG_NONNULL(3) SECP256K1_ARG_NONNULL(4) SECP256K1_ARG_NONNULL(5)
SECP256K1_ARG_NONNULL(8) SECP256K1_ARG_NONNULL(9);


/**
 *
 * @param ctx  - context object
 * @param gens - output
 * @param gen_seed - input seed
 * @return 1 - valid, 0 - invalid
 */
SECP256K1_API SECP256K1_WARN_UNUSED_RESULT int secp256k1_get_multi_generators(
        const secp256k1_context* ctx,
        secp256k1_generator *gens,
        const unsigned char *gen_seed,
        int m
) SECP256K1_ARG_NONNULL(1) SECP256K1_ARG_NONNULL(2) SECP256K1_ARG_NONNULL(3);

/**
 * Generate a multi-generator Pedersen commitment
 * @param ctx - context object
 * @param com - commitment
 * @param blind - secret key
 * @param digests - m x 32 bytes when each digest is 32 bytes (digests can be all zero if there is nothing to commit)
 * @param blind_gen - generator for the key
 * @param gens - generators for the values
 * @param m - number of maximum generators
 * @return
 */
SECP256K1_API SECP256K1_WARN_UNUSED_RESULT int secp256k1_get_multi_gen_commitment(
        secp256k1_pedersen_commitment *com,
        const unsigned char *blind,
        const unsigned char *digests,
        const secp256k1_generator *blind_gen,
        const secp256k1_generator *gens,
        int m
) SECP256K1_ARG_NONNULL(1) SECP256K1_ARG_NONNULL(2) SECP256K1_ARG_NONNULL(3) SECP256K1_ARG_NONNULL(4) SECP256K1_ARG_NONNULL(5);


SECP256K1_API void secp256k1_nopenena_context_create(
        const secp256k1_context* ctx,
        nopenena_context *nopenena_ctx,
        int n, // for short zerocom
        int maximum_p_size,
        int maximum_acc_size,
        const unsigned char *gen_seed
) SECP256K1_ARG_NONNULL(1) SECP256K1_ARG_NONNULL(2) SECP256K1_ARG_NONNULL(4);

void secp256k1_nopenena_context_clear(
        const secp256k1_context* ctx,
        nopenena_context *nctx
);

int secp256k1_create_account_sk(
        const secp256k1_context* ctx,
        nopenena_account_sk *acc_sk,
        unsigned char *r,
        unsigned char *k,
        unsigned char *gamma
) SECP256K1_ARG_NONNULL(1) SECP256K1_ARG_NONNULL(2) SECP256K1_ARG_NONNULL(3) SECP256K1_ARG_NONNULL(4) SECP256K1_ARG_NONNULL(5);


SECP256K1_API SECP256K1_WARN_UNUSED_RESULT int secp256k1_create_openout_sk(const secp256k1_context* ctx,
        nopenena_openout_sk *openout_sk,
        unsigned char *alpha,
        unsigned char *a,
        unsigned char *a_prime,
        unsigned char *b,
        unsigned char *kappa,
        unsigned char *rho
) SECP256K1_ARG_NONNULL(1) SECP256K1_ARG_NONNULL(2) SECP256K1_ARG_NONNULL(3) SECP256K1_ARG_NONNULL(4) SECP256K1_ARG_NONNULL(5) SECP256K1_ARG_NONNULL(6);

SECP256K1_API SECP256K1_WARN_UNUSED_RESULT int secp256k1_create_nopenena_account(
        const secp256k1_context* ctx,
        nopenena_context *nctx,
        nopenena_account *acc,
        nopenena_account_asset *asset,
        uint64_t value,
        nopenena_account_sk *acc_sk
) SECP256K1_ARG_NONNULL(1) SECP256K1_ARG_NONNULL(2) SECP256K1_ARG_NONNULL(3) SECP256K1_ARG_NONNULL(5);


SECP256K1_API int secp256k1_prove_nopenena_update(
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
        const unsigned char *nonce);


int secp256k1_verify_nopenena_update(
        const secp256k1_context* ctx,
        nopenena_context *nctx,
        nopenena_update_proof *proof,
        nopenena_account_asset *asset1,
        nopenena_account_asset *asset0,
        nopenena_account *acc,
        nopenena_openout *w);

SECP256K1_API  int secp256k1_open_nopenena_outputs(
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
        unsigned char *alpha_c_new);

void secp256k1_free_nopenena_outputs(
        nopenena_openout *out);

void secp256k1_free_nopenena_openout_proof(
        nopenena_openout_proof *proof);

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
        int zero_com_flag);

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
        int zero_com_flag);

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
        uint64_t withheld_c_new) ;

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
        secp256k1_pedersen_commitment *C_new);


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
);

int secp256k1_unlinked_logarithmic_zero_com_verify(
        const secp256k1_context* ctx,
        unsigned char *proof,
        int ring_size,
        int n,
        const unsigned char *challenge,
        const secp256k1_pedersen_commitment *coms,
        const secp256k1_generator *value_gen,
        const secp256k1_generator *blind_gen
);

#endif
