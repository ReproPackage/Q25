# Copyright 2025 <>
# Written by <>
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

library(tidyverse)
library(patchwork)
library(tikzDevice)
library(scales)
library(ggh4x)
library(ggpmisc)


INCH.PER.CM <- 1 / 2.54
WIDTH <- 18.1 * INCH.PER.CM

COL.WIDTH <- 8.85 * INCH.PER.CM
#COL.WIDTH <- 6.95 * INCH.PER.CM
COL.HEIGHT <- 6.2 * INCH.PER.CM

BASE.SIZE <- 9
FORMAT <- "tex"
theme_paper_base <- function() {
    return(theme_bw(base_size = BASE.SIZE) +
        theme(
            axis.title.x = element_text(size = BASE.SIZE),
            axis.title.y = element_text(size = BASE.SIZE),
            legend.title = element_text(size = BASE.SIZE),
            legend.position = "right",
            legend.box = "horizontal",
            legend.spacing.y = unit(-0.2, "cm"),
            plot.margin = unit(c(0.1, 0, 0, 0.1), "cm") #top, right, bottom, left
        ))
}

guide_paper_base <- function() {
    return(guides(
        shape = guide_legend(order = 1, nrow = 1, byrow = TRUE),
        col = guide_legend(order = 2, nrow = 10, byrow = TRUE, reverse = TRUE)
    ))
}

extract_monomial_density <- function(text) {
    matches <- gregexpr("\\d.\\d", text, perl= TRUE)
    extracted <- regmatches(text, matches)
    return(as.numeric(extracted))
}


OUT.PATH <- "./RPlots/"
ifelse(!dir.exists(file.path(OUT.PATH)),
        dir.create(file.path(OUT.PATH)),
        "")
POINT.ALPHA <- 0.6
POINT.SIZE <- 0.5
LINE.WIDTH <- 0.2

LFD.COLOURS <- c("black", "#E69F00", "#999999", "#009371", "#ed665a", "#beaed4", "#1f78b4", "#a44df0")
LFD.SHAPES <- c(15, 16, 17, 4, 5, 8, 9, 20)

###################################################################
####### Input data
###################################################################

df <- list.files(path="Experiments", pattern= "/*.csv", full.names=TRUE) %>% 
  lapply(read_csv) %>% 
  bind_rows

if ("qubo_percentile" %in% names(df))
{
    df <- mutate(df, qubo_percentile = ifelse(is.na(qubo_percentile), 1.0, qubo_percentile))
} else {
    df <- mutate(df, qubo_percentile=1.0)
}

relableMethod <- function(input){
    out = switch(input, KCMCQ={"Choi*"}, KDePFQ={"DeMorg."}, KT3CMCQ={"Choi"}, KT3DPFQ={"Dobry."})
    return(out)
}

relable_k <- function(orig) {
    return(paste("$k=", orig, "$", sep=""))
}

#Add SAT density:
df <- df %>% rowwise() %>%
    mutate(kSAT_maxUniqueClauses = choose(kSAT_actual_v, kSAT_k)) %>%
    mutate(kSAT_density = kSAT_c / kSAT_maxUniqueClauses) %>%
    filter(kSAT_density <= 1.0) %>%
    filter(kSAT_c %in% c(288, 236, 184, 132, 93, 54, 15)) %>%
    filter(kSAT_k %in% c(3,7,10,15,20)) %>% 
    rowwise() %>%
    mutate(Path=relableMethod(Path), kSAT_k=relable_k(kSAT_k))

df$kSAT_kf = factor(df$kSAT_k, levels=c("$k=3$","$k=7$","$k=10$","$k=15$","$k=20$")) 
    
dfP1 <- filter(df, qubo_percentile==1.0)

# ###################################################################
# ####### Plot ksat variables vs qubo variables
# ###################################################################

plt.variables <- ggplot(dfP1, aes(x=kSAT_actual_v, y=qubo_variables, colour=as.factor(kSAT_c))) +
    facet_grid(cols=vars(kSAT_kf), rows=vars(Path))+
    geom_point(size=POINT.SIZE, alpha=POINT.ALPHA) +
    geom_line(linewidth = LINE.WIDTH) +
    theme_paper_base() +
    guide_paper_base() +
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x))) +
    annotation_logticks(sides = "l") +
    scale_colour_manual("\\# Clauses", values=LFD.COLOURS) +
    xlab("\\# Variables in k-SAT") +
    ylab("\\# Variables in QUBO [log]")
ggsave(plot=plt.variables, filename="Variables.pdf", path=OUT.PATH, width=2*COL.WIDTH, height=COL.HEIGHT)


# ###################################################################
# ####### Plot variables vs runtime
# ###################################################################

plt.VarsvsTime <- ggplot(dfP1, aes(x=kSAT_actual_v, y=time_totalPath, colour=as.factor(kSAT_c))) +
    facet_grid(cols=vars(kSAT_kf), rows=vars(Path))+
    geom_point(size=POINT.SIZE, alpha=POINT.ALPHA) +
    geom_line(linewidth = LINE.WIDTH) +
    theme_paper_base() +
    guide_paper_base() +
    scale_y_log10(breaks = c(0.01, 1, 100), labels = trans_format("log10", math_format(10^.x))) +
    annotation_logticks(sides = "l") +
    scale_colour_manual("\\# Clauses", values=LFD.COLOURS) +
    xlab("\\# Variables in k-SAT") +
    ylab("Time for path transformations [s, log]")
ggsave(plot=plt.VarsvsTime, filename="VarsvsTime.pdf", path=OUT.PATH, width=2*COL.WIDTH, height=COL.HEIGHT)


# ###################################################################
# ####### Plot density vs density
# ###################################################################

plt.Density <- ggplot(dfP1, aes(x=kSAT_density, y=qubo_density_deg2, colour=as.factor(kSAT_c))) +
    facet_grid(cols=vars(kSAT_kf), rows=vars(Path))+
    geom_point(size=POINT.SIZE, alpha=POINT.ALPHA) +
    theme_paper_base() +
    guide_paper_base() +
    scale_colour_manual("\\# Clauses", values=LFD.COLOURS) +
    scale_shape_manual("\\# Variables", values=LFD.SHAPES) +
    xlab("k-SAT density") +
    ylab("Qubo density degree 2")
ggsave(plot=plt.Density, filename="Density.pdf", path=OUT.PATH, width=2*COL.WIDTH, height=COL.HEIGHT)

# ###################################################################
# ####### Plot variables vs monomials in qubo
# ###################################################################

plt.variablesMonomials <- ggplot(dfP1, aes(x=kSAT_actual_v, y=qubo_monomials, colour=as.factor(kSAT_c))) +
    facet_grid(cols=vars(kSAT_kf), rows=vars(Path))+
    geom_point(size=POINT.SIZE, alpha=POINT.ALPHA) +
    geom_line(linewidth = LINE.WIDTH) +
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x))) +
    annotation_logticks(sides = "l") +
    theme_paper_base() +
    guide_paper_base() +
    scale_colour_manual("\\# Clauses", values=LFD.COLOURS) +
    xlab("\\# Variables in k-SAT") +
    ylab("\\# Monomials in QUBO [log]")
ggsave(plot=plt.variablesMonomials, filename="VariablesMonomials.pdf", path=OUT.PATH, width=2*COL.WIDTH, height=COL.HEIGHT)


# ###################################################################
# ####### Plot variables vs variables in qubo for KDePFQ
# ###################################################################

plt.variablesPercentile <- ggplot(filter(df, Path == "DeMorg." & qubo_percentile %in% c(0.01, 0.1, 1.0)), aes(x=kSAT_actual_v, y=qubo_variables, colour=as.factor(kSAT_c))) +
    facet_grid(cols=vars(kSAT_kf), rows=vars(qubo_percentile))+
    geom_point(size=POINT.SIZE, alpha=POINT.ALPHA) +
    geom_line(linewidth = LINE.WIDTH) +
    theme_paper_base() +
    guide_paper_base() +
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x))) +
    annotation_logticks(sides = "l") +
    scale_colour_manual("\\# Clauses", values=LFD.COLOURS) +
    xlab("\\# Variables in k-SAT") +
    ylab("\\# Variables in QUBO [log]")
ggsave(plot=plt.variablesPercentile, filename="VariablesPercentile.pdf", path=OUT.PATH, width=2*COL.WIDTH, height=COL.HEIGHT)
