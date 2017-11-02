## initialize packrat (install same version of packages that was used by the
# authors and saves them to a local library)

## load needed libraries
library(ggplot2)
theme_set(theme_bw())
library(magrittr)
library(tidyr)
library(purrr)
library(dplyr)
library(survival)
library(mgcv)
devtools::install_github("adibender/pammtools@v0.0.3.1")
library(pammtools)
library(grid)
library(gridExtra)

## color pallets
Set1    <- RColorBrewer::brewer.pal(9, "Set1")
Greens  <- RColorBrewer::brewer.pal(9, "Greens")
Purples <- RColorBrewer::brewer.pal(9, "Purples")


#### Motivational figure
hweibull <- function(t, alpha=1.5, tau=2) {
  (alpha/tau) * (t/tau)^(alpha -1)
}

tdf <- data.frame(time = seq(0, 4, by=0.05)) %>%
  mutate(hweib = hweibull(time))

set.seed(16042017)
sdata <- data.frame(time=rweibull(500, 1.5, 2)) %>%
  mutate(event=1)

gg.base <- ggplot(tdf, aes(x=time, y=hweib)) +
  geom_line() +
  ylab(expression(lambda(t))) +
  xlab("t")



kappa5 <- seq(0, 4, by=0.8)
ped5 <- split_data(Surv(time, event)~., data=sdata, cut=kappa5, id="id")
pam5 <- gam(ped_status ~ interval, data=ped5, family=poisson(), offset=offset)
pred5 <- ped5 %>% int_info() %>% add_hazard(pam5)

pdf("weibullHazard.pdf", width=9, height=3)
grid.draw(
    cbind(
        ggplotGrob(gg.base),
        ggplotGrob(gg.base +  geom_vline(xintercept = kappa5, lty=2) +
          scale_x_continuous(
                breaks = kappa5,
            labels=c(
              expression(kappa[0]==0),
              expression(kappa[1]==0.8),
              expression(kappa[2]==1.6),
              expression(kappa[3]==2.4),
              expression(kappa[4]==3.2),
              expression(kappa[5]==4)))),
        ggplotGrob(gg.base + geom_vline(xintercept = kappa5, lty=2) +
          geom_step(data=plot_df(pred5), aes(x=time, y=hazard)) +
          scale_x_continuous(
                breaks = kappa5,
            labels=c(
              expression(kappa[0]==0),
              expression(kappa[1]==0.8),
              expression(kappa[2]==1.6),
              expression(kappa[3]==2.4),
              expression(kappa[4]==3.2),
              expression(kappa[5]==4)))),
        size = "last"))
dev.off()



#### Initialize veterans data
data("veteran", package="survival")

# we limit the data to event times < 400 (not many events afterwards), for
# convenience
veteran %<>%
  mutate(
    prior = 1*(prior!=0),
    trt   = 1*(trt==2)) %>%
  filter(time < 400)
head(veteran)


## creation of example table
library(xtable)
set.seed(2017)
index <- sample(1:nrow(veteran))
vet_tex <- xtable(veteran[index[1:6],], digits = rep(0,9))
print(vet_tex, include.rownames=FALSE)


#### Baseline

## basic nelson-aalen estimate
base.cox <- basehaz(coxph(Surv(time, status)~1, data=veteran))

## data trafo for pem/pam
vet_ped <- split_data(Surv(time, status)~., data = veteran, id   = "id")
head(vet_ped)

## veteran data in ped output for paper
vet_ped_tex <- xtable::xtable(
  vet_ped[1:6,-c(8:11)],
  digits=c(rep(0,5),1,rep(0,3)))
print(vet_ped_tex,include.rownames=FALSE)

## basic survival models PEM and PAM
pem_obj <- gam(ped_status ~ interval, data = vet_ped, offset = offset, family = poisson())
pam_obj <- gam(ped_status ~ s(tend), data  = vet_ped, offset = offset, family = poisson())

## extract time/interval information for plotting
int_df <- int_info(vet_ped)
int_df <- int_df %>%
  add_cumu_hazard(pam_obj) %>%
  left_join(base.cox, by=c("tend"="time"))
int_df$pemhaz <- predict(pem_obj, newdata=int_df, type="response")
int_df %<>% mutate(pemch = cumsum(pemhaz * intlen))

## Figure comparing Nelson-Aalen, PEM and PAM estimates
pdf("compar_cumhaz.pdf",width = 7,height = 4)
ggplot(int_df, aes(x = tend)) +
  geom_ribbon(aes(ymin = cumu_lower, ymax = cumu_upper), alpha = 0.2) +
  geom_step(aes(y = hazard,    col = "Nelson-Aalen")) +
  geom_line(aes(y = pemch,     col = "PEM")) +
  geom_line(aes(y = cumu_hazard, col = "PAM")) +
  scale_color_manual(
    name   = "Method",
    values = c("PEM" = Set1[2],"PAM" = 1,"Nelson-Aalen" = Set1[1])) +
  theme(legend.position = "bottom") +
  ylab(expression(hat(Lambda)(t))) + xlab("t") +
  ggtitle("Comparison of cumulative hazards using Nelson-Aalen vs. PEM vs. PAM")
dev.off()

## difference between pem and nelson-allen (evaluated at interval end points)
all.equal(int_df$hazard, int_df$pemch)


#### Model with time-constant effects of time-constant covariates

## Cox PH
cph <- coxph(Surv(time, status) ~ trt + pspline(karno, df = 0), data = veteran)
## PAM
pam <- gam(ped_status ~ s(tend) + trt + s(karno, bs = "ps"), data = vet_ped,
  family = poisson(), offset = offset)

# create data set with varying karno values (from min to max in 100 steps)
# and add term contribution of karno from PAM and Cox PH models
karno_df <- vet_ped %>%
  make_newdata(expand = "karno", n = 100) %>%
  add_term(pam, term = "karno") %>%
  mutate(
    cox = predict(cph, newdata = ., type = "terms")[, "pspline(karno, df = 0)"] -
      predict(cph, newdata=data.frame(trt=mean(veteran$trt), karno=mean(veteran$karno))))

# Figure comparing smooth effect estimates of Cox PH and PAM models
pdf("compar_smooth.pdf", width=4, height=2)
ggplot(karno_df, aes(x=karno, ymin=low, ymax=high)) +
  geom_line(aes(y=fit, col="PAM")) +
  geom_ribbon(alpha=0.2) +
  geom_line(aes(y=cox, col="Cox PH"))+
  scale_colour_manual(name="Method",values=c("Cox PH"=Set1[1],"PAM"=1)) +
  xlab(expression(x[plain(karno)])) + ylab(expression(hat(f)(x[plain(karno)])))
dev.off()


#### Stratified Proportional Hazards Model
table(veteran$celltype)

## stratified Cox PH
cph  <- coxph(Surv(time, status) ~ strata(celltype), data=veteran)
strata_df <- basehaz(cph)

## prepare plot object comparing stratified Cox PH and stratified PAM
gg.bz <- ggplot(strata_df, aes(x=time, y=hazard, group=strata)) +
  geom_step(aes(col="stratified Cox")) +
  ylab(expression(hat(Lambda)(t))) + xlab("t") + facet_wrap(~strata)

## stratified PAM (note to include main effects of celltype)
# using first difference penalties here
pam <- gam(ped_status ~ celltype + s(tend, by=celltype),
  data=vet_ped, family=poisson(), offset=offset)

## prepare data set for plotting
pinf <- vet_ped %>%
  group_by(celltype) %>%
  ped_info() %>%
  add_cumu_hazard(pam) %>%
  rename(strata=celltype) %>%
  filter(!(strata == "adeno" & tend > 200))


pdf("stratifiedPH.pdf", width=6, height=3)
gg.bz +
    geom_line(data=pinf, aes(x=tend, y=cumu_hazard, group=strata, col="stratified PAM")) +
    geom_ribbon(data=pinf, aes(x=tend, y=cumu_hazard, ymin=cumu_lower, ymax=cumu_upper), alpha=0.2) +
    facet_wrap(~strata, nrow=1) +
    scale_color_manual(name="Method", values=c(Set1[1], 1)) +
    theme(legend.position="bottom")
dev.off()


#### Linear, smoothly time-varying effects
head(veteran)
vfit <- coxph(Surv(time, status) ~ trt + prior + karno + tt(karno),
    data = veteran,
    tt = function(x, t, ...) x * log(t + 20))
coef(vfit)
ttcoef <- round(coef(vfit), 3)[3:4]
t <- seq(0, 400, by = 10)

# data transformation
vet_ped <- vet_ped %>%
  mutate(logt20 = log((tstart)+(tstart-tend)/2 + 20))
head(vet_ped) %>% select(interval, ped_status, trt, karno, age, prior, logt20)

# fit PAM with log-transform
pam <- gam(ped_status ~ s(tend) + trt + prior + karno + karno:logt20,
    data = vet_ped, offset = offset, family = poisson())
coef(pam)[1:5]

# fit penalized PAM (no need to specify main effect for karno here)
pam2 <- gam(ped_status ~ s(tend) + trt + prior + s(tend, by = karno),
    data = vet_ped, offset = offset, family = poisson())
term_df <- vet_ped %>%
  ped_info() %>%
  add_term(pam2, term = "karno") %>%
  mutate_at(c("fit", "low", "high"), ~ ./karno) %>%
  mutate(
    cox.fit = coef(vfit)["karno"] + coef(vfit)["tt(karno)"]*log(tend + 20),
    pam.fit = coef(pam)["karno"]  + coef(pam)["karno:logt20"]*log(tend + 20))

pdf("linearSmoothTVkarno.pdf", width=5, height=2)
ggplot(term_df, aes(x = tend, y = fit)) +
    geom_step(aes(col="PAM with penalized spline")) +
    geom_stepribbon(aes(ymin = low, ymax = high), alpha = 0.2) +
    geom_line(aes(y = cox.fit, col = "Cox with log-transform")) +
    geom_step(aes(y = pam.fit, col = "PAM with log-transform")) +
    scale_color_manual(name="Method", values = c(Set1[1:2], "black")) +
    xlab("t") + ylab(expression(hat(f)(t)))
dev.off()



#### Non-linear, non-linearly time-varying effects

# PAM with tensor product smooth
pam3 <- gam(ped_status ~ trt + prior + s(age) + te(tend, karno),
  data=vet_ped, family=poisson(), offset=offset)

## heat map/contour plot
te_gg <- gg_tensor(pam3) +
  geom_vline(xintercept = c(1, 51, 200), lty = 3) +
  geom_hline(yintercept = c(40, 75, 95), lty = 3) +
  scale_fill_gradient2(
    name = expression(hat(f)(list(x[plain(karno)],t))),
    low  = "steelblue", high = "firebrick2") +
  geom_contour(col="grey30") +
  xlab("t") + ylab(expression(x[plain(karno)])) +
  theme(
    axis.text        = element_text(size = rel(1.2)),
    axis.title       = element_text(size = rel(1.3)),
    legend.position  = "bottom",
    legend.title     = element_text(size = rel(1.3)),
    legend.text      = element_text(size = rel(1.2)),
    strip.background = element_blank(),
    strip.text.x     = element_blank())


## plot f(karno, t) for specific slices
karno_df <- combine_df(
  int_info(vet_ped),
  select(sample_info(vet_ped), -karno),
  data.frame(karno = c(40, 75, 95)))

karno_df <- karno_df %>% add_term(pam3, term="karno")

karno_gg <- ggplot(karno_df, aes(x=tend, y=fit)) +
  geom_step(aes(col=factor(karno)), lwd=1.1) +
  geom_stepribbon(aes(ymin = low, ymax = high, fill=factor(karno)), alpha=.2) +
  scale_color_manual(
    name   = expression(x[plain(karno)]),
    values = Greens[c(4, 7, 9)]) +
  scale_fill_manual(
    name   = expression(x[plain(karno)]),
    values = Greens[c(4, 7, 9)]) +
  ylab(expression(hat(f)(list(x[plain(karno)],t)))) +
  xlab("t")+  coord_cartesian(ylim= c(-4, 3)) +
  theme(
    axis.text       = element_text(size = rel(1.2)),
    axis.title      = element_text(size = rel(1.3)),
    legend.position = "bottom",
    legend.title    = element_text(size = rel(1.3)),
    legend.text     = element_text(size = rel(1.2)))

time_df <- vet_ped %>%
  make_newdata(tend=c(1, 51, 200), karno=seq(20, 100, by=5)) %>%
  add_term(pam3, term="karno")

time_gg <- ggplot(time_df, aes(x=karno)) +
  geom_line(aes(y=fit, col=factor(tend)), lwd=1.1) +
  geom_ribbon(aes(ymin = low, ymax = high, fill=factor(tend)), alpha=.2) +
  scale_color_manual(name="t", values=Purples[c(4, 6, 8)]) +
  scale_fill_manual(name="t", values=Purples[c(4, 6, 8)]) +
  ylab(expression(hat(f)(list(x[plain(karno)],t)))) +
  xlab(expression(x[plain(karno)])) + coord_cartesian(ylim= c(-4, 3)) +
  theme(
    axis.text       = element_text(size = rel(1.2)),
    axis.title      = element_text(size = rel(1.3)),
    legend.position = "bottom",
    legend.title    = element_text(size = rel(1.3)),
    legend.text     = element_text(size = rel(1.2)))

library(grid)
library(gridExtra)
pdf("karnoTEgarphs.pdf", width=9, height=4)
grid.arrange(
 te_gg,
 karno_gg,
 time_gg, nrow=1)
dev.off()


####  recidivism data: Time-dependent covariates

# raw data
# http://socserv.mcmaster.ca/jfox/Books/Companion/scripts/appendix-cox.R
prison <- read.table("http://math.unm.edu/~james/Rossi.txt", header = TRUE) %>%
  mutate(subject=row_number())

# transform into long format (equivalent to time-discrete models b/c the offset
# is not needed here (log(1)=0)
prison_long <- prison %>%
  gather(calendar.week, employed, emp1:emp52) %>%
  filter(!is.na(employed)) %>% # employed unequal to NA only for intervals under risk
  group_by(subject) %>%
  mutate(
    start = row_number()-1,
    stop = row_number(),
    arrest = ifelse(stop==last(stop) & arrest==1, 1, 0),
    offset = log(stop - start)) %>%
  select(subject, start, stop, offset, arrest, employed, fin:educ) %>%
  arrange(subject, stop)
prison_long %<>% mutate(employed.lag1 = lag(employed, default=0)) %>%
  slice(-1) %>% # this makes data set equivalent to analysis by John Fox
  ungroup()

psamp <- filter(prison_long, subject %in% 1:2) %>%
  select(subject, start, stop, arrest, employed, fin, age, mar) %>% print(n=100)
xtab_long <- xtable::xtable(psamp, digits=rep(0, 9))
print(xtab_long, include.rownames=FALSE)
# fit pem with smoth baseline
# note that gam call for both models completely equivalent except for offset,
# which, however, is always 0 in this case, because interval length is always
# exactly one week and event time are discrete
pam <- gam(arrest ~ s(stop) + fin + s(age, bs="ps") + race + wexp + mar + paro +
  s(prio, bs="ps") + employed.lag1,
  data=prison_long, family=poisson(), offset=offset)
summary(pam)
cph <- coxph(Surv(start, stop, arrest)~ fin + pspline(age) + race + wexp + mar +
  paro + pspline(prio) + employed.lag1,
  data=prison_long)
summary(cph)

## prep coefficient plot
# hideous hack from here (rotates error bar legend):
# http://stackoverflow.com/questions/40229263/how-to-rotate-legend-symbols-in-ggplot2
GeomPointrange$draw_key <-  function (data, params, size)     {

         draw_key_vpath <- function (data, params, size) {
           # only need to change the x&y coords so that the line is horizontal
           # originally, the vertical line was `0.5, 0.1, 0.5, 0.9`
              segmentsGrob(0.1, 0.5, 0.9, 0.5,
              gp = gpar(
                col = alpha(data$colour, data$alpha),
                lwd = data$size * .pt,
                lty = data$linetype,
                lineend = "butt"),
              arrow = params$arrow)
              }

    grobTree(draw_key_vpath(data, params, size),
             draw_key_point(transform(data, size = data$size * 4), params))
}

# -> very similar coefficient estimates  (and baseline estimates)
exp(cbind(coef(pam)[2:7], coef(cph)[c(1, 14:17, 30)]))

all_eff <- purrr::map_df(
  list(
    tidy_fixed(pam),
    tidy_fixed(cph)[-c(2:3, 8:9), ]),
  bind_rows, .id="Model") %>%
  mutate(Model = factor(Model, levels=2:1, labels=c("Cox-PH", "PAM")))


## coefficient plot (fixed effects) recidivism data
gg.coef <- ggplot(all_eff, aes(x=variable, y=coef, ymin=lower, ymax=upper)) +
  geom_hline(yintercept = 0, lty=3) +
  geom_pointrange(aes(col=factor(Model), shape=factor(Model)),
    position=position_dodge(width=0.5)) +
  scale_colour_manual(name="Method", values=c(1, Set1[1]), limits=rev(levels(all_eff$Model))) +
  scale_shape_manual(name="Method",values=c(19, 15), limits=rev(levels(all_eff$Model))) +
  ylab(expression(hat(beta)%+-% 1.96 %.% SE)) +
  xlab("") +
  coord_flip(ylim=range(-1.5, 1)) +
  theme(legend.position="bottom")

pinf <- prison_long %>%
  make_newdata(expand="age") %>%
  add_term(pam, term="age") %>%
  mutate(cphfit = predict(object=cph, ., type="terms")[,2])


gg.age <- ggplot(pinf, aes(x=age, y=fit)) +
  geom_line(aes(col="PAM")) +
  geom_ribbon(aes(ymin=low, ymax=high), alpha=0.2) +
  geom_line(aes(y=cphfit, col="Cox-PH")) +
  scale_colour_manual(name="Method", values=c(Set1[1], 1), breaks=c("Cox-PH", "PAM")) +
  xlab("Age") + ylab(expression(hat(f)(x))) +
  theme(legend.position="none")
## prio

pinf.prio <- prison_long %>%
  make_newdata(expand="prio") %>%
  add_term(pam, term="prio") %>%
  mutate(cphfit = predict(object=cph, ., type="terms")[,"pspline(prio)"])

gg.prio <- ggplot(pinf.prio, aes(x=prio, y=fit)) +
  geom_line(aes(col="PAM")) +
  geom_ribbon(aes(ymin=low, ymax=high), alpha=0.2) +
  geom_line(aes(y=cphfit, col="Cox-PH")) +
  scale_colour_manual(name="Method", values=c(Set1[1], 1)) +
  xlab("Number of prior convictions") +
  ylab(expression(hat(f)(x))) +
  theme(legend.position="none")


pdf("prisonAllEffects.pdf", width=6, height=4)
grid.arrange(
  gg.coef + theme(legend.position="bottom"),
  gg.age,
  gg.prio,
  layout_matrix=matrix(c(1, 1, 2, 3), ncol=2))
dev.off()
