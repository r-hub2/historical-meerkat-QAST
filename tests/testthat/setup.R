# Variancias (jackknife/bootstrap) no processo principal durante os testes:
# o covr so registra a cobertura do processo principal, e os refits rodam
# nos filhos do fork quando ha mais de um worker. O caminho com fork segue
# testado explicitamente em test-jackknife-se.R (cores = 2).
options(QAST.cores = 1L)
