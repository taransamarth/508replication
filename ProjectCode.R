library(tidyverse)
library(haven)
library(igraph)
library(MASS)
library(gridExtra)
library(statnet)
library(assortnet)

# Network construction

state_court_cites <- read_dta("HinkleNelson_Files/data.for.analysis.dta") %>%
  rename(weight = citeCount)

state_net <- graph_from_data_frame(state_court_cites %>% filter(weight > 0), directed = T)

E(state_net)$arrow.size <- .05

set.seed(2022)

GGally::ggnet2(state_net, edge.size = E(state_net)$weight / 10, node.label = names(V(state_net)))

ggsave("FinalProjectFigs/StateCiteNetwork.pdf", width = 5, height = 5, units = "in")

# Histogram for vertex attributes

vertex_attrs <- read_dta("HinkleNelson_Files/StateLevelVariables_1.dta") %>%
  left_join(read_dta("HinkleNelson_Files/cfScores.dta")) %>%
  left_join(read_dta("HinkleNelson_Files/berryData.dta") %>% dplyr::select(-c(year))) %>% filter(state != "DC")

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

# Network descriptive statistics

### Centralities

cent_stats <- data.frame(
  in_deg = degree(state_net, mode = "in"),
  out_deg = degree(state_net, mode = "out"),
  in_strength = strength(state_net, mode = "in"),
  out_strength = strength(state_net, mode = "out"),
  eigen = evcent(state_net, directed = TRUE)$vector
) %>% rownames_to_column("state")

top_10 <- data.frame(
  in_deg = cent_stats$state[(sort(cent_stats$in_deg, index.return = T, decreasing = T))$ix],
  out_deg = cent_stats$state[(sort(cent_stats$out_deg, index.return = T, decreasing = T))$ix],
  in_strength = cent_stats$state[(sort(cent_stats$in_strength, index.return = T, decreasing = T))$ix],
  out_strength = cent_stats$state[(sort(cent_stats$out_strength, index.return = T, decreasing = T))$ix],
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

