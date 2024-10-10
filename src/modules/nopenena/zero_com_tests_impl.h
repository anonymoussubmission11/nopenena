/**********************************************************************
 * Copyright (c) 2024 Jayamine Alupotha                                 *
 * Distributed under the MIT software license, see the accompanying   *
 * file COPYING2 or http://www.opensource.org/licenses/mit-license.php.*
 **********************************************************************/

#ifndef RAHAS_SECP256K1_TESTS_IMPL_H
#define RAHAS_SECP256K1_TESTS_IMPL_H

#include "secp256k1_nopenena.h"
#include "openssl/rand.h"
#include "math.h"

void run_identity_tests(void) {
    secp256k1_pedersen_commitment com;
    secp256k1_pedersen_commitment *commitments;
    secp256k1_unlinked_identity_pf proof1;
    unsigned char blind[32];
    unsigned char nonce[32];
    secp256k1_hash digest;
    int i, t;
    int index, index_of_com, index_of_gen, ring_size;
    double linear_time_prove = 0;
    double log_time_prove = 0;
    double linear_time_verify = 0;
    double log_time_verify = 0;
    time_t start;
    time_t end;
    int m = 5;

    secp256k1_generator *generators = malloc(m * sizeof(secp256k1_generator));
    unsigned char *digests = (unsigned char*) malloc(m * 32);
    unsigned char **blinds = (unsigned char**) malloc(10 * sizeof(unsigned char*));
    for (i = 0; i < 10; i++)
        blinds[i] = (unsigned char*) malloc(32);

    /* Commitment creation */
    secp256k1_rand256(blind);
    secp256k1_rand256(nonce);
    secp256k1_rand256(digest.data);
    CHECK(secp256k1_identity_create(ctx, &com, blind, digest, &secp256k1_generator_const_h, &secp256k1_generator_const_g));

    /* unlinked proving */
    ring_size = 10;

    commitments = (secp256k1_pedersen_commitment *) malloc(ring_size * sizeof(secp256k1_pedersen_commitment));
    for (i = 0; i < ring_size; i++) {
        secp256k1_rand256(blinds[i]);
        digest.data[0] = i;
        CHECK(secp256k1_identity_create(ctx, &commitments[i], blinds[i], digest, &secp256k1_generator_const_h,
                                        &secp256k1_generator_const_g));
    }

    for (index = 0; index < 10; index++) {
        digest.data[0] = index;
        start = clock();
        CHECK(secp256k1_unlinked_logarithmic_identity_prove(ctx, &proof1, index, blinds[index],
                                                            digest.data, nonce, 10, 4, commitments,
                                                            &secp256k1_generator_const_h, &secp256k1_generator_const_g) == 1);
        end = clock();
        log_time_prove += ((double) (end - start)) / CLOCKS_PER_SEC;

        start = clock();
        CHECK(secp256k1_unlinked_logarithmic_identity_verify(ctx, &proof1, digest.data, nonce, commitments, 10, 4,
                                                             &secp256k1_generator_const_h, &secp256k1_generator_const_g) == 1);
        end = clock();
        log_time_verify += ((double) (end - start)) / CLOCKS_PER_SEC;

        free(proof1.data);
    }
    printf("linear prove: %f\n", linear_time_prove);
    printf("linear verify: %f\n", linear_time_verify);
    printf("log prove: %f\n", log_time_prove);
    printf("log verify: %f\n", log_time_verify);

    index_of_com = 5;
    secp256k1_rand256(nonce);
    for (index_of_gen = 0; index_of_gen < m; index_of_gen++) {
        digest.data[0] = index_of_com;
        CHECK(secp256k1_unlinked_logarithmic_identity_prove(ctx, &proof1, index_of_com, blinds[index_of_com],
                                                            digest.data, nonce, 10, 4, commitments,
                                                            &secp256k1_generator_const_h, &secp256k1_generator_const_g) == 1);

        CHECK(secp256k1_unlinked_logarithmic_identity_verify(ctx, &proof1, digest.data, nonce, commitments, 10, 4,
                                                             &secp256k1_generator_const_h, &secp256k1_generator_const_g) == 1);

        free(proof1.data);
    }


    /* Multi Generator Proofs */
    m = 3;
    memset(digests, 0, m*32);
    for (t = 0; t < m; t++) {
        digests[t*32] = t + 1;
    }

    CHECK(secp256k1_get_multi_generators(ctx, generators, nonce, m));
    for (i = 0; i < ring_size; i++) {
        CHECK(secp256k1_get_multi_gen_commitment(&commitments[i],
                                                 blinds[i],
                                                 digests,
                                                 &secp256k1_generator_const_g,
                                                 generators, m));
    }

    free(digests);
    for (i = 0; i < ring_size; i++) {
        free(blinds[i]);
    }
    free(blinds);
    free(generators);
    free(commitments);

}



#endif /* RAHAS_SECP256K1_TESTS_IMPL_H */