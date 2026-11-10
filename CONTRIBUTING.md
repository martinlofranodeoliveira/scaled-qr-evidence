# Contributing

Contributions to `scaled-qr-evidence` are welcome within its documented scope.

1. Run offline unit and LAPACK oracle tests before submitting changes:
   ```sh
   /bin/sh scripts/test-offline.sh
   /bin/sh scripts/demo-offline.sh --assert
   ```
2. Keep the production solver free of external BLAS/LAPACK runtime linkage.
3. General contact: [Lofrano Tecnologia](https://lofranotecnologia.com.br).
