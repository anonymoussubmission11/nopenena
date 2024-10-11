# Nopenena Untraceable Payments

This projects modifies Libsecp256k1 library with functions for Nopenena.

## Installation and Testing

```
$ ./autogen.sh
$ ./configure
$ make
$ ./tests
$ sudo make install  # optional
```

## Nopenena Untraceable Payment Examples

We give typical Nopenena payments, Nopenena payments with an excrow mechansim, and Nopenena split payments in module [link](https://github.com/anonymoussubmission11/nopenena/tree/main/src/modules/nopenena).

The following code shows the benchmark code for typical Nopenane payments,

```C
void bench_nopenena_api_t(int zero_com_flag, int n) {
    int P[] = {2, 2, 2, 2, 2, 2};
    int A[] = {32, 8, 12, 16, 20, 24};
    for (int test = 0; test < 6; test++) {
        int p_size = P[test]; // power of two only
        int acc_size = A[test];
        int all_combo_tests = 2;
        nopenena_context nctx;
        nopenena_account acc[acc_size];
        nopenena_account_asset asset0[acc_size];
        nopenena_account_asset asset1[acc_size];
        nopenena_update_proof uproof[acc_size];
        nopenena_account_sk acc_sk[acc_size];
        nopenena_openout_sk openout_sk[acc_size];
        nopenena_account_asset *asset0_ptr[acc_size];
        nopenena_account_asset *asset1_ptr[acc_size];
        nopenena_update_proof *uproof_ptr[acc_size];
        nopenena_openout_sk *openout_sk_ptr[acc_size];
        nopenena_openout_proof proof;
        nopenena_openout w;
        uint64_t values[acc_size];
        int flags[acc_size];
        int j_index[p_size];
        uint64_t values_new[acc_size];
        uint64_t values_ptr[p_size];
        uint64_t values_new_ptr[p_size];
        nopenena_account_sk *acc_sk_ptr[p_size];
        nopenena_balance_proof pi_balance;
        unsigned char r[32];
        unsigned char k[32];
        unsigned char rho[32];
        unsigned char gamma[32];
        unsigned char alpha[32];
        unsigned char a[32];
        unsigned char a_prime[32];
        unsigned char b[32];
        unsigned char kappa[32];
        unsigned char gen_seed[32];
        unsigned char nonce[32];
        unsigned char nonce1[32];
        uint64_t withheld_c = 0;
        uint64_t withheld_c_new = 0;
        int i;

        memset(gen_seed, 0, 32);
        memset(nonce, 1, 32);

        secp256k1_nopenena_context_create(ctx, &nctx, n, p_size, acc_size, gen_seed);

        secp256k1_rand256(rho);

        // create accounts
        for (i = 0; i < acc_size; i++) {
            flags[i] = 0;
            values[i] = 40;
            values_new[i] = values[i];

            // Create secrets
            secp256k1_rand256(r);
            secp256k1_rand256(k);
            secp256k1_rand256(gamma);
            CHECK(secp256k1_create_account_sk(ctx, &acc_sk[i], r, k, gamma) == 1);
            secp256k1_rand256(alpha);
            secp256k1_rand256(a);
            secp256k1_rand256(a_prime);
            secp256k1_rand256(b);
            secp256k1_rand256(kappa);
            CHECK(secp256k1_create_openout_sk(ctx, &openout_sk[i], alpha, a, a_prime, b, kappa, rho) == 1);

            // Create accounts
            CHECK(secp256k1_create_nopenena_account(ctx, &nctx, &acc[i], &asset0[i], values[i], &acc_sk[i]) == 1);
        }

        for (i = 0; i < acc_size; i++) {
            asset0_ptr[i] = &asset0[i];
            asset1_ptr[i] = &asset1[i];
            openout_sk_ptr[i] = &openout_sk[i];
            uproof_ptr[i] = &uproof[i];
        }


        int l = 0;

        // All combo ======================================================
        for (; l < all_combo_tests; l++) {

            if (p_size != 3) {
                // real participants
                for (l = 0; l < p_size / 2; l++) {
                    j_index[l] = l;
                    values_new[j_index[l]] = values[j_index[l]] + 1;
                    flags[j_index[l]] = 1;
                }
                for (l = p_size / 2; l < p_size; l++) {
                    j_index[l] = l;
                    values_new[j_index[l]] = values[j_index[l]] - 1;
                    flags[j_index[l]] = 1;
                }
            } else {
                for (l = 0; l < p_size; l++) {
                    j_index[l] = l;
                    flags[j_index[l]] = 1;
                }
                values_new[j_index[0]] = values[j_index[0]] + 1;
                values_new[j_index[1]] = values[j_index[1]] + 1;
                values_new[j_index[2]] = values[j_index[2]] - 2;
            }

            secp256k1_rand256(r);
            secp256k1_rand256(nonce);

            time_t start = clock();
            for (i = 0; i < acc_size; i++) {
                CHECK((flags[i] == 0 && values_new[i] - values[i] == 0) ||
                      (flags[i] == 1 && values_new[i] - values[i] != 0));

                // Update the randomness
                CHECK(secp256k1_create_account_sk_r(ctx, &acc_sk[i], r) == 1);
                CHECK(secp256k1_update_nopenena_account(ctx, &nctx, asset1_ptr[i], asset0_ptr[i], &acc[i], values_new[i],
                                                      values[i],
                                                      &acc_sk[i]) == 1);
            }

            // create pointers
            for (i = 0; i < p_size; i++) {
                values_ptr[i] = values[j_index[i]];
                values_new_ptr[i] = values_new[j_index[i]];
                acc_sk_ptr[i] = &acc_sk[j_index[i]];
            }

            CHECK(secp256k1_open_nopenena_outputs(ctx, &nctx, &w,
                                                asset1_ptr,
                                                values_new_ptr,
                                                values_ptr,
                                                acc_sk_ptr,
                                                openout_sk_ptr,
                                                j_index, p_size, acc_size, nonce,
                                                withheld_c, NULL, withheld_c_new, NULL) == 1);

            for (i = 0; i < acc_size; i++) {
                CHECK(secp256k1_prove_nopenena_update(ctx, &nctx, &uproof[i],
                                                    asset1_ptr[i], asset0_ptr[i], &acc[i],
                                                    values_new[i], values[i], &acc_sk[i], &openout_sk[i],
                                                    &w, flags[i], i, nonce) == 1);
            }

            secp256k1_rand256(nonce);

            CHECK(secp256k1_prove_nopenena_outputs(ctx, &nctx, &proof, &w,
                                                 asset1_ptr,
                                                 values_new_ptr,
                                                 values_ptr,
                                                 acc_sk_ptr,
                                                 openout_sk_ptr,
                                                 j_index, p_size, acc_size, uproof_ptr, nonce, NOPENI_ALL,
                                                 zero_com_flag) == 1);

            secp256k1_rand256(nonce);
            secp256k1_rand256(nonce1);
            CHECK(secp256k1_prove_nopenena_balance(ctx, &nctx, &pi_balance, 0, 0, acc, asset0_ptr, asset1_ptr,
                                                 values, values_new,
                                                 acc_size, r, nonce, nonce1, NULL, NULL, 0, NULL, NULL, 0) == 1);
            time_t end = clock();
            double prove_all = ((double) (end - start)) / CLOCKS_PER_SEC;

            double verify_all = 0;
            for (int trial = 0; trial < 3; trial++) {
                start = clock();

                for (i = 0; i < acc_size; i++) {
                    CHECK(secp256k1_verify_nopenena_update(ctx, &nctx, &uproof[i], \
                                             asset1_ptr[i], asset0_ptr[i], &acc[i], &w) == 1);
                }
                CHECK(secp256k1_verify_nopenena_outputs(ctx, &nctx, &proof, &w,
                                                      asset1_ptr, p_size, acc_size, uproof_ptr, NOPENI_ALL,
                                                      zero_com_flag) == 1);
                CHECK(secp256k1_verify_nopenena_balance(ctx, &nctx, &pi_balance, 0, 0, acc, asset0_ptr, asset1_ptr,
                                                      acc_size, NULL, NULL) == 1);

                end = clock();
                verify_all += ((((double) (end - start)) / CLOCKS_PER_SEC)*1000);
            }

            long int txsize_without_pi =
                    w.pi_len + acc_size * (4 * 32 + 2 * 33) + p_size * (3 * 33 + 2 * 32) + 3 * 33 +
                    32 + 32;
            printf("%d, %d, %d, %d, %d, %f, %f, %d ,%ld, %ld, %ld\n", acc_size, p_size,
                   secp256k1_nopenena_get_total_combo(p_size, acc_size), nctx.rctx.n, nctx.rctx.m,
                   prove_all, verify_all/3, proof.pi_zero_len, w.pi_len, txsize_without_pi, proof.pi_zero_len + txsize_without_pi);

            secp256k1_free_nopenena_outputs(&w);
            secp256k1_free_nopenena_openout_proof(&proof);
            for (l = 0; l < p_size; l++) {
                flags[j_index[l]] = 0;
                values[j_index[l]] = values_new[j_index[l]];
            }
            for (i = 0; i < acc_size; i++) { // swap them
                if (l % 2 == 0) {
                    asset0_ptr[i] = &asset1[i];
                    asset1_ptr[i] = &asset0[i];
                } else {
                    asset0_ptr[i] = &asset0[i];
                    asset1_ptr[i] = &asset1[i];
                }
            }
        }

        secp256k1_nopenena_context_clear(ctx, &nctx);
    }
}
```

