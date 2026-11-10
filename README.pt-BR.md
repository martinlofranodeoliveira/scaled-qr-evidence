# Scaled QR Evidence (PT-BR)

Versão de desenvolvimento `0.1.0` — biblioteca e CLI de mínimos quadrados densos em Fortran 2008 (`gfortran` / `fpm`) com fatoração QR Householder original com escalonamento e pivoteamento de colunas, até 3 passos de refinamento iterativo e diagnósticos numéricos verificados contra o oráculo SVD nativo `DGELSD` do LAPACK.

Desenvolvido por [Lofrano Tecnologia](https://lofranotecnologia.com.br). [English version](README.md).

## Visão Geral

- **Domínio**: matrizes `real64` finitas com $1 \le n \le 32$, $n \le m \le 128$, $1 \le k_{\text{rhs}} \le 4$, magnitudes não nulas em $[2^{-100}, 2^{100}]$ e limiar de triagem $2^{-k}$ ($k \in [10, 40]$).
- **Resultados tipados**: `COMPUTED`, `INVALID_INPUT`, `NUMERIC_UNRESOLVED`, `UNKNOWN_BUDGET`.
- **Oráculo de testes**: compara soluções em problemas bem-condicionados contra `DGELSD` (diferença relativa $\le 10^{-8}$) sem vincular o LAPACK ao binário de produção.

## Como Executar

```sh
/bin/sh scripts/test-offline.sh
/bin/sh scripts/demo-offline.sh --assert
/bin/sh scripts/demo-offline.sh --assert --json
```

## Limitações e Exclusões de Escopo

- Triagem de pivôs é uma política numérica de ponto flutuante, não um certificado algébrico formal.
- Não resolve sistemas subdeterminados ($m < n$) nem problemas esparsos.

## Licença

MIT — Copyright (c) 2026 Lofrano Tecnologia. Consulte [LICENSE](LICENSE).
