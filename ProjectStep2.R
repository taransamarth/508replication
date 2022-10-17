library(tidyverse)
library(haven)
library(igraph)
library(MASS)
library(gridExtra)

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
  left_join(read_dta("HinkleNelson_Files/berryData.dta") %>% dplyr::select(-c(year)))

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
