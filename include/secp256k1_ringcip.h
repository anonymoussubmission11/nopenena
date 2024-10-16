//
// Created by jayamine on 10/24/23.
//

#ifndef SECP256K_NOPENI_SECP256K1_RINGCIP_H
#define SECP256K_NOPENI_SECP256K1_RINGCIP_H

typedef struct ringcip_context_struct {
    int L; // range bit
    int n; // N = n^m when N is the ring size and n is the base
    int m;
    int N;
    secp256k1_ge *multigen;
    secp256k1_generator geng;
    secp256k1_generator genh;
    secp256k1_generator genmu;
} ringcip_context;

typedef struct cint_public_struct {
    secp256k1_ge c; // could be uninitialized
    uint8_t buf[33]; // could be uninitialized
} cint_pt;

typedef struct cint_secret_struct {
    int64_t v;  // could be uninitialized
    secp256k1_scalar val;  // could be uninitialized
    secp256k1_scalar key; // could be uninitialized
} cint_st;

typedef struct zero_com_proof_struct {
    secp256k1_ge A; // could be uninitialized
    uint8_t bufA[32]; // could be uninitialized
    secp256k1_ge B; // could be uninitialized
    uint8_t bufB[32]; // could be uninitialized
    secp256k1_ge C; // could be uninitialized
    uint8_t bufC[32]; // could be uninitialized
    secp256k1_ge D; // could be uninitialized
    uint8_t bufD[32]; // could be uninitialized
} gb_zero_proof_t;


/**
 * Create the context object
 * @param ctx - main context
 * @param L -  range bit
 * @param n -
 * @param m -
 * @param gen_seed - generate seed
 * @return ring cip context
 */
SECP256K1_API  ringcip_context secp256k1_ringcip_context_create(
        const secp256k1_context* ctx,
        int L, int n, int m, uint8_t *gen_seed,
        secp256k1_generator *blind_gen)
SECP256K1_ARG_NONNULL(1);

/**
 * Clear the ring cip context
 * @param ctx - main context
 * @param rctx - ring cip context
 */
SECP256K1_API  void secp256k1_ringcip_context_clear(ringcip_context *rctx)
SECP256K1_ARG_NONNULL(1) SECP256K1_ARG_NONNULL(2);

/**
 * create a confidential integer
 * @param ctx
 * @param rctx - ring cip context
 * @param c  - output
 * @param v - hidden integer
 * @param key - blind factor
 * @return 1- success, 0 - failure
 */
int secp256k1_create_secret_cint(const secp256k1_context* ctx,
                          ringcip_context *rctx,
                          cint_st *c,
                          int64_t v,
                          uint8_t *key
 ) SECP256K1_ARG_NONNULL(1) SECP256K1_ARG_NONNULL(2) SECP256K1_ARG_NONNULL(3) SECP256K1_ARG_NONNULL(5);

/**
 * create a confidential integer
 * @param ctx
 * @param rctx - ring cip context
 * @param c  - output
 * @param csk  - input secrets
 * @return 1- success, 0 - failure
 */
int secp256k1_create_cint(const secp256k1_context* ctx,
                                 ringcip_context *rctx,
                                 cint_pt *c,
                                 cint_st *csk
) SECP256K1_ARG_NONNULL(1) SECP256K1_ARG_NONNULL(2) SECP256K1_ARG_NONNULL(3) SECP256K1_ARG_NONNULL(5);


/**
 * Get serialized commitment
 * @param ctx - context
 * @param buf - output bytes [33]
 * @param c - confidential integer
 * @return 1- successful, 0 - failure
 */
int secp256k1_serialize_cint(
        const secp256k1_context* ctx, uint8_t *buf, cint_pt *c)
SECP256K1_ARG_NONNULL(1) SECP256K1_ARG_NONNULL(2) SECP256K1_ARG_NONNULL(3);

/**
 * Get serialized commitment
 * @param ctx - context
 * @param c - confidential integer
 * @param buf - output bytes [33]
 * @return 1- successful, 0 - failure
 */
int secp256k1_parse_cint(const secp256k1_context* ctx, cint_pt *c, uint8_t *buf)
SECP256K1_ARG_NONNULL(1) SECP256K1_ARG_NONNULL(2) SECP256K1_ARG_NONNULL(3);

/**
 * Commitment of multiple generators
 * @param ctx - context
 * @param rctx - ring cip context
 * @param com - output
 * @param vals - values
 * @param key - random key
 * @return 1- successful, 0 - failure
 */
int secp256k1_create_multival_com(const secp256k1_context* ctx,
                                  const ringcip_context *rctx,
                                  secp256k1_ge *com,
                                  secp256k1_scalar *vals,
                                  secp256k1_scalar *key)
SECP256K1_ARG_NONNULL(1) SECP256K1_ARG_NONNULL(2) SECP256K1_ARG_NONNULL(3) SECP256K1_ARG_NONNULL(4) SECP256K1_ARG_NONNULL(5);

/**
 * return the byte size of the proof
 * @param rctx
 * @return
 */
int secp256k1_zero_mcom_get_size(const ringcip_context *rctx);

/**
 * create a proof according to
``Short Accountable Ring Signatures Based on DDH'' by
 Jonathan Bootle, Andrea Cerulli, Pyrros Chaidos, Essam Ghadafi, Jens Groth & Christophe Petit
 * @param ctx - context object
 * @param rctx - ring cip context
 * @param proof - output (please allocate memory for the proof)
 * @param Cs -  commitment set
 * @param index - index of the zero value commitment
 * @param key - secret key of the zero value commitment
 * @return 1- successful, 0 - failure
 */
int secp256k1_create_zero_mcom_proof(const secp256k1_context* ctx,
                                     const ringcip_context *rctx,
                                     uint8_t *proof,
                                     cint_pt *Cs,
                                     int index,
                                     secp256k1_scalar *key)
SECP256K1_ARG_NONNULL(1) SECP256K1_ARG_NONNULL(2) SECP256K1_ARG_NONNULL(3) SECP256K1_ARG_NONNULL(4) SECP256K1_ARG_NONNULL(6);

/**
 * verify a proof according to
``Short Accountable Ring Signatures Based on DDH'' by
 Jonathan Bootle, Andrea Cerulli, Pyrros Chaidos, Essam Ghadafi, Jens Groth & Christophe Petit
 * @param ctx - context object
 * @param rctx - ring cip context
 * @param proof - input
 * @param Cs commitment set
 * @return 1- successful, 0 - failure
 */
int secp256k1_verify_zero_mcom_proof(const secp256k1_context* ctx,
                                     const ringcip_context *rctx,
                                     uint8_t *proof,
                                     cint_pt *Cs)
SECP256K1_ARG_NONNULL(1) SECP256K1_ARG_NONNULL(2) SECP256K1_ARG_NONNULL(3) SECP256K1_ARG_NONNULL(4);

#endif //SECP256K_NOPENI_SECP256K1_RINGCIP_H


