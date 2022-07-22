probsDip <- expression(rbind(
  dnbinom(txmin:txmax, tkcov/tbias, mu = tkcov * 1),
  dnbinom(txmin:txmax, tkcov/tbias, mu = tkcov * 2)
))

probsTri <- expression(rbind(
  dnbinom(txmin:txmax, tkcov/tbias, mu = tkcov * 1),
  dnbinom(txmin:txmax, tkcov/tbias, mu = tkcov * 2),
  dnbinom(txmin:txmax, tkcov/tbias, mu = tkcov * 3)
))

probsTet <- expression(rbind(
  dnbinom(txmin:txmax, tkcov/tbias, mu = tkcov * 1),
  dnbinom(txmin:txmax, tkcov/tbias, mu = tkcov * 2),
  dnbinom(txmin:txmax, tkcov/tbias, mu = tkcov * 3),
  dnbinom(txmin:txmax, tkcov/tbias, mu = tkcov * 4)
))



factorDip <- expression(
  c(
    2*tth/(1+tth),
    1/(1+tth)
  )
)
factorAut <- expression(
  c(4*tth/(3+tth),
    6*tth/(6+5*tth+tth^2),
    8*tth/(6+11*tth+6*tth^2+tth^3),
    6/(6+11*tth+6*tth^2+tth^3))
)
factorAll <- expression(
  c(
    4*exp(-3*tth*tdiverg)*tth*(
      -2*exp(tdiverg*(tth-2)) +
        exp(3*tdiverg*tth) * (2+tth)^2 * (3 + tth) -
        2*exp(2*tth*tdiverg)*(3+4*tth+tth^2)
    ),
    2*exp(-1/2*tdiverg*(4+5*tth))*(
      6*exp(tdiverg*tth/2)*tth+
        exp(1/2*tdiverg*(4+5*tth)) * (2+tth)^2 * (3 + tth) +
        2*exp(2*tdiverg+3/2*tdiverg*tth)*
        (-6-8*tth+tth^2+tth^3)
    ),
    8*exp(-tdiverg*(3+4*tth))*tth*(
      -exp(tdiverg+2*tdiverg*tth)+
        exp(3+tdiverg*(1+tth))*(3+tth)
    ),
    2*exp(-2*tdiverg*(1+tth))*(
      tth + 2*exp(tdiverg*(2+tth))*(3+tth)
    )
  )/ (1+tth) / (2+tth)^2 / (3+tth)

)

factorTse <- expression(
  pal * c(
    4*exp(-3*tth*tdiverg)*tth*(
      -2*exp(tdiverg*(tth-2)) +
        exp(3*tdiverg*tth) * (2+tth)^2 * (3 + tth) -
        2*exp(2*tth*tdiverg)*(3+4*tth+tth^2)
    ),
    2*exp(-1/2*tdiverg*(4+5*tth))*(
      6*exp(tdiverg*tth/2)*tth+
        exp(1/2*tdiverg*(4+5*tth)) * (2+tth)^2 * (3 + tth) +
        2*exp(2*tdiverg+3/2*tdiverg*tth)*
        (-6-8*tth+tth^2+tth^3)
    ),
    8*exp(-tdiverg*(3+4*tth))*tth*(
      -exp(tdiverg+2*tdiverg*tth)+
        exp(3+tdiverg*(1+tth))*(3+tth)
    ),
    2*exp(-2*tdiverg*(1+tth))*(
      tth + 2*exp(tdiverg*(2+tth))*(3+tth)
    )
  )/ (1+tth) / (2+tth)^2 / (3+tth) +
    (1-pal) *   c(4*tth/(3+tth),
                  6*tth/(6+5*tth+tth^2),
                  8*tth/(6+11*tth+6*tth^2+tth^3),
                  6/(6+11*tth+6*tth^2+tth^3))
)

factorTraaa <- expression(
  c(
    3*tth / (2 + tth),
    3*tth / (2 + 3*tth + tth^2),
    2 / (2 + 3*tth + tth^2)
  )
)

factorTraab <- expression(
  c(
    exp(-tdiverg*tth)*(-2-4*tth+exp(tdiverg*tth)*(2+7*tth+3*tth^2)) / (2+3*tth+tth^2),

    (1 + ((2*exp(-tdiverg*tth) * (tth-1))/(2+tth)))/ (1 + tth),

    2*exp(-tdiverg*tth) / (2 + 3*tth + tth^2)
  )
)
