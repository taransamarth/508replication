library(tidyverse)
library(haven)
library(igraph)
library(MASS)
library(gridExtra)
library(statnet)
library(assortnet)

### Helpers
# Function for AMEN compatibility with broom and huxtable
tidy.ame <- function(model, term_names = colnames(model$BETA)) {
  output <- tibble(data.frame(
    term = term_names,
    estimate = colMeans.mcmc.list(model$BETA),
    std.error = apply(model$BETA, 2, sd)) %>% mutate(
      statistic = estimate/std.error,
      p.value = 2*(1-pnorm(abs(statistic)))))
  
  return(output)
}

### Replication

# Network construction

state_court_cites <- read_dta("HinkleNelson_Files/data.for.analysis.dta") %>%
  rename(weight = citeCount)

state_net <- graph_from_data_frame(state_court_cites %>% filter(weight > 0), directed = T)

E(state_net)$arrow.size <- .05

set.seed(2022)

GGally::ggnet2(state_net, edge.size = E(state_net)$weight / 20, arrow.gap = .03,
               arrow.size = 3,
               node.label = names(V(state_net)))

ggsave("FinalProjectFigs/StateCiteNetwork.pdf", width = 5, height = 5, units = "in")

# Histogram for vertex attributes

vertex_attrs <- read_dta("HinkleNelson_Files/StateLevelVariables_1.dta") %>%
  left_join(read_dta("HinkleNelson_Files/cfScores.dta")) %>%
  left_join(read_dta("HinkleNelson_Files/berryData.dta") %>% dplyr::select(-c(year))) %>% filter(state != "DC")

geo_edge_attrs <- read_dta("HinkleNelson_Files/ivDyads.dta") %>% 
  left_join(read_dta("HinkleNelson_Files/migrationData.dta")) %>%
  filter(citingCourt != "DC" & citedCourt != "DC")

sq_hist <-ggplot(vertex_attrs) + geom_histogram(aes(x = SquireIndex)) + xlab("Squire Professionalism Index")
ggsave("FinalProjectFigs/SquireHistogram.pdf")

case_hist <- ggplot(vertex_attrs) + geom_histogram(aes(x = Caseload)) + xlab("Yearly Caseload")
ggsave("FinalProjectFigs/CaseloadHistogram.pdf")

pop_hist <- ggplot(vertex_attrs) + geom_histogram(aes(x = Pop2010)) + xlab("2010 State Population")
ggsave("FinalProjectFigs/PopHistogram,pdf")

elect_hist <- ggplot(vertex_attrs) + geom_histogram(aes(x = Elected)) + xlab("Elected Court")
ggsave("FinalProjectFigs/ElectedHistogram.pdf")

tort_hist <- ggplot(vertex_attrs) + geom_histogram(aes(x = TortScore)) + xlab("Tort Score")
ggsave("FinalProjectFigs/TortHistogram.pdf")

legcap_hist <- ggplot(vertex_attrs) + geom_histogram(aes(x = LegCap)) + xlab("Legal Capital")
ggsave("FinalProjectFigs/LegalCapitalHistogram.pdf")

cfscore_hist <- ggplot(vertex_attrs) + geom_histogram(aes(x = cfscore)) + xlab("Median Judge Ideology (Bonica cfScore)")
ggsave("FinalProjectFigs/MedianJudgeHistogram.pdf")

berry_hist <- ggplot(vertex_attrs) + geom_histogram(aes(x = citi6013)) + xlab("Citizen Ideology (Berry)")
ggsave("FinalProjectFigs/CitizenIdeologyHistogram.pdf")

gridarr <- grid.arrange(sq_hist, case_hist, 
                        pop_hist, elect_hist, 
                        tort_hist, legcap_hist, 
                        cfscore_hist, berry_hist)

ggsave(plot = gridarr, filename = "FinalProjectFigs/Histograms.pdf")

# Regression replication

rep_model <- glm.nb(weight ~ ctIDdistance + berryCitDif + bornInCitedState + distanceS + diffPopMil + as.factor(geobucket) +
  sameMethod + citedSquireProf + citedLegCapS + citedPopMil + citedElected2 +
  citingElected2 + citedLA + citingLA + s2sCites2 + citedTortScore, data = state_court_cites)

# Regression table

stargazer::stargazer(rep_model, style = "apsr", covariate.labels = c(
  "ID Distance: Courts",
  "ID Distance: Citizens",
  "Cultural linkage",
  "Distance", "Pop. difference",
  "Same West region only",
  "Same federal circuit only",
  "Contiguous only",
  "Same West and circuit",
  "Same West and contiguous",
  "Same circuit and contiguous",
  "Same West, circuit, and contiguous",
  "Same selection method",
  "Cited ct: legal professionalism",
  "Cited ct: Legal capital",
  "Cited ct: Population",
  "Cited ct: Elected",
  "Citing ct: Elected",
  "Cited ct: Louisiana",
  "Citing ct: Louisiana",
  "Total cites difference",
  "Cited ct: Tort score"
))

# Adding vertex attributes to network

V(state_net)$Caseload <- vertex_attrs$Caseload
V(state_net)$SquireIndex <- vertex_attrs$SquireIndex
V(state_net)$Pop2010 <- vertex_attrs$Pop2010
V(state_net)$Method <- as.factor(vertex_attrs$Method)
V(state_net)$Elected <- as.factor(vertex_attrs$Elected)
V(state_net)$LegCap <- vertex_attrs$LegCap
V(state_net)$Circuit <- as.factor(vertex_attrs$Circuit)
V(state_net)$cfscore <- vertex_attrs$cfscore
V(state_net)$citi6013 <- vertex_attrs$citi6013
V(state_net)$TortScore <- vertex_attrs$TortScore

contig_matrix <- as.matrix(geo_edge_attrs %>% 
                             dplyr::select(-c(distance, bornInCitedState)) %>% 
                             spread(citingCourt, contiguous))
rownames(contig_matrix) <- colnames(contig_matrix)[-1]
contig_matrix <- contig_matrix[,-1]
mode(contig_matrix) <- "numeric"
contig_matrix[is.na(contig_matrix)] <- 0

dist_matrix <- as.matrix(geo_edge_attrs %>% 
                             dplyr::select(-c(contiguous, bornInCitedState)) %>% 
                             spread(citingCourt, distance))
rownames(dist_matrix) <- colnames(dist_matrix)[-1]
dist_matrix <- dist_matrix[,-1]
mode(dist_matrix) <- "numeric"
dist_matrix[is.na(dist_matrix)] <- 0

born_matrix <- as.matrix(geo_edge_attrs %>% 
                           dplyr::select(-c(contiguous, distance)) %>% 
                           spread(citingCourt, bornInCitedState))
rownames(born_matrix) <- colnames(born_matrix)[-1]
born_matrix <- born_matrix[,-1]
mode(born_matrix) <- "numeric"
born_matrix[is.na(born_matrix)] <- 0

# Network descriptive statistics

### Centralities

cent_stats <- data.frame(
  in_deg = igraph::degree(state_net, mode = "in"),
  out_deg = igraph::degree(state_net, mode = "out"),
  in_strength = igraph::strength(state_net, mode = "in"),
  out_strength = igraph::strength(state_net, mode = "out"),
  eigen = igraph::evcent(state_net, directed = TRUE)$vector
) %>% rownames_to_column("state")

top_10 <- data.frame(
  in_deg = cent_stats$state[(sort(cent_stats$in_deg, index.return = T, decreasing = T))$ix],
  out_deg = cent_stats$state[(sort(cent_stats$out_deg, index.return = T, decreasing = T))$ix],
  in_strength = cent_stats$state[(sort(cent_stats$in_strength, index.return = T, decreasing = T))$ix],
  out_strength = cent_stats$state[(sort(cent_stats$out_strength, index.return = T, decreasing = T))$ix],
  in_strength_w = cent_stats$state[(sort(cent_stats$in_strength/V(state_net)$LegCap, index.return = T, decreasing = T))$ix],
  out_strength_w = cent_stats$state[(sort(cent_stats$out_strength/V(state_net)$Caseload, index.return = T, decreasing = T))$ix],
  eigen = cent_stats$state[(sort(cent_stats$eigen, index.return = T, decreasing = T))$ix]
)

knitr::kable(head(top_10, 10), format = "latex")

cor(V(state_net)$LegCap, cent_stats$out_strength)
cor(V(state_net)$LegCap, cent_stats$out_deg)
cor(V(state_net)$LegCap, cent_stats$in_strength)
cor(V(state_net)$LegCap, cent_stats$in_deg)

### Reciprocity

el2 <- reshape2::melt(adj_mat)
el2$y <- reshape2::melt(t(adj_mat))$value

cor(el2$value, el2$y, use = "complete.obs")

### Assortativity

in_out <- assortativity(state_net, strength(state_net, mode = "in"), 
              strength(state_net, mode = "out"))

out_in <- assortativity(state_net, strength(state_net, mode = "out"), 
              strength(state_net, mode = "in"))

in_in <- assortativity(state_net, strength(state_net, mode = "in"), 
              strength(state_net, mode = "in"))

out_out <- assortativity(state_net, strength(state_net, mode = "out"), 
              strength(state_net, mode = "out"))

assort_cf <- assortment.continuous(as_adjacency_matrix(state_net, attr = "weight"), 
                    vertex_values = V(state_net)$cfscore, SE = T)

assort_berry <- assortment.continuous(as_adjacency_matrix(state_net, attr = "weight"), 
                      vertex_values = V(state_net)$citi6013, SE = T)

assort_circ <- assortment.discrete(as_adjacency_matrix(state_net, attr = "weight"), 
                    types = as.factor(V(state_net)$Circuit))

knitr::kable(data.frame(
  types = c("Strength (In, Out)", "Strength (Out, In)", "Strength (In, In)", "Strength (Out, Out)",
            "Court Ideology", "Citizen Ideology", "Circuit"),
  coefs = round(c(in_out, out_in, in_in, out_out, assort_cf$r, assort_berry$r, assort_circ$r),3)
), format = "latex")


# Modeling

## Prepare covariates

diff.CitizenIdeo <- as.matrix(dist(V(state_net)$citi6013, upper = TRUE))
diff.CourtIdeo <- as.matrix(dist(V(state_net)$cfscore, upper = TRUE))
diff.Circuit <- as.matrix(dist(V(state_net)$Circuit, upper = TRUE)) + 1
diff.Circuit[diff.Circuit != 1] <- 0
diff.Method <- as.matrix(dist(V(state_net)$Method, upper = TRUE)) + 1
diff.Method[diff.Method != 1] <- 0
diff.Population <- as.matrix(log(dist(V(state_net)$Pop2010, upper = TRUE)))
total.Cites <- matrix(rowSums(as.matrix(state_net[,])), nrow = 52, ncol = 52) - as.matrix(state_net[,])

covar.Squire <- as.matrix((V(state_net)$SquireIndex))
covar.CourtIdeo <- as.matrix((V(state_net)$citi6013))
covar.CitizenIdeo <- as.matrix((V(state_net)$cfscore))
covar.TortScore <- as.matrix((V(state_net)$TortScore))
covar.LegCap <- log(as.matrix((V(state_net)$LegCap)))
covar.Caseload <- log(as.matrix((V(state_net)$Caseload)))
covar.Pop2010 <- log(as.matrix((V(state_net)$Pop2010)))
covar.Elected <- as.matrix(V(state_net)$Elected) - 1

Xdyad <- abind(Circuit = diff.Circuit, 
               Citizen = diff.CitizenIdeo, 
               CourtIdeo = diff.CourtIdeo,
               Contiguous = contig_matrix,
               Distance = dist_matrix,
               Born = born_matrix,
               Method = diff.Method,
               TotalCites = total.Cites,
               along = 3)

Xs <- data.frame(Caseload = covar.Caseload,
                 Elected = covar.Elected,
                 TortScore = covar.TortScore,
                 LegCap = covar.LegCap)

Xr <- data.frame(Squire = covar.Squire,
                 LegCap = covar.LegCap, 
                 Population = covar.Pop2010,
                 Elected = covar.Elected,
                 Caseload = covar.Caseload,
                 TortScore = covar.TortScore)

## Run AMEN model
amen_model2 <- ame(as.matrix(state_net[,]), 
                  Xdyad = Xdyad, 
                  Xrow = as.matrix(Xs),
                  Xcol = as.matrix(Xr),
                  R = 2, symmetric = FALSE, 
                  family = "nrm", rvar = TRUE, cvar = TRUE, print = TRUE, nscan = 100000,
                  burn = 1000, odens = 100, seed = 2022)

amen_terms <- c("Intercept", "Sender.Caseload", "Sender.Elected", "Sender.TortScore", 
                "Sender.LegCap", "Receiver.Squire", "Receiver.LegCap", "Receiver.Pop2010",
                "Receiver.Elected", "Receiver.Caseload", "Receiver.TortScore",
                "Dyad.Circuit", "Dyad.Citizen", "Dyad.CourtIdeo", "Dyad.Contiguous",
                "Dyad.Distance", "Dyad.Born", "Dyad.Method", "Dyad.TotalCites")

## Run ERGM models
net_state <- intergraph::asNetwork(state_net)

set.vertex.attribute(net_state, "LegCap", log(get.vertex.attribute(net_state, "LegCap")))
set.vertex.attribute(net_state, "Pop2010", log(get.vertex.attribute(net_state, "Pop2010")))

job::job({ergm_ocovar <- ergm(net_state ~ sum + nonzero + 
                            absdiff("cfscore") + absdiff("citi6013") + edgecov(born_matrix, attrname = "Born") + 
                            edgecov(dist_matrix, attrname = "Distance") + nodematch("Circuit") + 
                            edgecov(contig_matrix, attrname = "Contiguous") + nodematch("Method") + nodeicov("SquireIndex") + 
                            nodeicov("LegCap") + nodeicov("Pop2010") + nodeicov("Elected") + 
                            nodeocov("Elected") + edgecov(total.Cites, attrname = "Total Cites") + nodeicov("TortScore") + 
                            nodeocov("LegCap") + nodeocov("TortScore") + nodeocov("Caseload") + 
                            nodeicov("Caseload") + mutual(form = "min") + CMP + transitiveweights("min", "max", "min") +
                              nodeocovar(transform = "sqrt"), 
                          response = "weight", reference = ~Poisson, verbose = 0,
                          control = control.ergm(seed = 2022, MCMC.samplesize = 4000))})

job::job({ergm_icovar <- ergm(net_state ~ sum + nonzero + 
                                absdiff("cfscore") + absdiff("citi6013") + edgecov(born_matrix, attrname = "Born") + 
                                edgecov(dist_matrix, attrname = "Distance") + nodematch("Circuit") + 
                                edgecov(contig_matrix, attrname = "Contiguous") + nodematch("Method") + nodeicov("SquireIndex") + 
                                nodeicov("LegCap") + nodeicov("Pop2010") + nodeicov("Elected") + 
                                nodeocov("Elected") + edgecov(total.Cites, attrname = "Total Cites") + nodeicov("TortScore") + 
                                nodeocov("LegCap") + nodeocov("TortScore") + nodeocov("Caseload") + 
                                nodeicov("Caseload") + mutual(form = "min") + CMP + transitiveweights("min", "max", "min") +
                                nodeicovar(transform = "sqrt"), 
                              response = "weight", reference = ~Poisson, verbose = 0,
                              control = control.ergm(seed = 2022, MCMC.samplesize = 4000))})

ergm_terms <- c("sum", "nonzero", "Dyad.CourtIdeo", "Dyad.Citizen", "Dyad.Born",
                "Dyad.Distance", "Dyad.Circuit", "Dyad.Contiguous", "Dyad.Method",
                "Receiver.Squire", "Receiver.LegCap", "Receiver.Pop2010", "Receiver.Elected",
                "Sender.Elected", "Dyad.TotalCites", "Receiver.TortScore", "Sender.LegCap",
                "Sender.TortScore", "Sender.Caseload", "Receiver.Caseload", "mutual.min",
                "CMP", "transitiveweights.min.max.min")

## Prepare results
results <- huxreg(list("AMEN" = tidy_override(amen_model2, term = amen_terms),
                       "ERGM (ocovar)" = tidy_override(ergm_ocovar, term = c(ergm_terms, "sqrt(nodeocovar)"), extend = TRUE), 
                       "ERGM (icovar)" = tidy_override(ergm_icovar, term = c(ergm_terms, "sqrt(nodeicovar)"), extend = TRUE)), 
                  statistics = c("AIC", N = "nobs"))

quick_latex(results, file = "models.tex", open = TRUE)

## Goodness-of-fit

ergm_sims <- simulate(ergm_ocovar, nsim=1000, seed = 2022, response = "weight")

### Compute estimates and drop simulation object
ergm_stats <- data.frame(do.call(rbind, 
                                 lapply(ergm_sims, function(x){
                                   gofstats(as.matrix(x, attrname = "weight"))})))

rm(ergm_sims)

ergm_sims <- simulate(ergm_icovar, nsim=1000, seed = 2022, response = "weight")

### Compute estimates and drop simulation object
ergm_stats2 <- data.frame(do.call(rbind, 
                                  lapply(ergm_sims, function(x){
                                    gofstats(as.matrix(x, attrname = "weight"))})))

rm(ergm_sims)

### Retrieve dependency estimates from AMEN fits
amen_stats <- data.frame(amen_model2$GOF)

### Retrieve actual dependency values from network
net_stats <- rownames_to_column(data.frame("actual" = gofstats(as.sociomatrix(net_state, attrname = "weight"))))

### Plot dependency values to compare across models
gof_df <- rbind(pivot_longer(amen_stats, colnames(amen_stats)) %>% 
                  mutate(model = "AMEN \n(nrm)"), 
                pivot_longer(ergm_stats, colnames(ergm_stats)) %>% 
                  mutate(model = "ERGM \n(sender)"), 
                pivot_longer(ergm_stats2, colnames(ergm_stats2)) %>% 
                  mutate(model = "ERGM \n(receiver)")) %>% 
  left_join(net_stats, by = c("name" = "rowname")) %>% 
  group_by(name, model) 

gof_df %>% 
  summarise(mean = mean(value), high = quantile(value, .95), 
            low = quantile(value, .05), high90 = quantile(value, .90), 
            low10 = quantile(value, .1), actual = unique(actual)) %>%
  ggplot(aes(x = as.factor(model))) + geom_hline(aes(yintercept = actual)) +
  geom_pointrange(aes(ymin = low, y = mean, ymax = high), color = "#5A5A5A", alpha = .7) +
  geom_pointrange(aes(ymin = low10, y = mean, ymax = high90)) +
  facet_wrap(~name, scales = "free", nrow = 1) + xlab("Model") + ylab("Parameter Estimate")

### Save plot
ggsave("ModelComparison.pdf", units = "in", width = 15, height = 5)

### Save dataset
write.csv(gof_df, "gof_stats.csv")

### Convergence plots

mcmc.diagnostics(ergm_icovar)
mcmc.diagnostics(ergm_ocovar)
plot(amen_model2)



