all_freq<-read.table("NAC055_insertion_all_freq_CAN.txt")
colnames(all_freq) <- c("pop", "freq", "n")

all_freq <- all_freq[order(all_freq$pop),]

climate_data <- climate_data[climate_data$ID %in% all_freq$pop,]

all_freq$AI_july <-climate_data$AI_Monthly.awi_pm_sr_07/10000
all_freq$AI_August <-climate_data$AI_Monthly.awi_pm_sr_08/10000


all_freq <- as.data.frame(all_freq)
2/(2*sum(all_freq$n[all_freq$pop %in% c("ES01", "ES02", "ES03", "ES06", "ES10")]))

sum(all_freq$n[all_freq$pop %in% c("ES01", "ES02", "ES03", "ES06", "ES10")])
df <- all_freq

# Calculate counts
df$ALT <- round(df$freq * 2 * df$n)
df$REF <- 2 * df$n - df$ALT

# Logistic regression using allele counts
model <- glm(cbind(ALT, REF) ~ AI_july, data = df, family = binomial)

# Summary of the model
summary(model)


library(ggplot2)

# Add predicted values to data frame
df$predicted <- predict(model, type = "response")

ggplot(df, aes(x = AI_july, y = freq)) +
  geom_point(size = 2) +
  geom_line(aes(y = predicted), color = "blue", size = 1) +
  labs(x = "AI (July)", y = "Allele Frequency", title = "Logistic Regression Fit") +
  theme_minimal()

model <- glm(cbind(ALT, REF) ~ AI_August, data = df, family = binomial)

summary(model)
# Add predicted values to data frame
df$predicted <- predict(model, type = "response")

ggplot(df, aes(x = AI_August, y = freq)) +
  geom_point(size = 2) +
  geom_line(aes(y = predicted), color = "blue", size = 1) +
  labs(x = "AI (July)", y = "Allele Frequency", title = "Logistic Regression Fit") +
  theme_minimal()




###


df <- all_freq[!all_freq$pop%in%c("ES02", "ES10", "ES15"),]
df$ALT <- round(df$freq * 2 * df$n)
df$REF <- 2 * df$n - df$ALT


# Logistic regression using allele counts
model <- glm(cbind(ALT, REF) ~ AI_july, data = df, family = binomial)

# Summary of the model
summary(model)


library(ggplot2)

# Add predicted values to data frame
df$predicted <- predict(model, type = "response")

ggplot(df, aes(x = AI_july, y = freq)) +
  geom_point(size = 2) +
  geom_line(aes(y = predicted), color = "blue", size = 1) +
  labs(x = "AI (July)", y = "Allele Frequency", title = "Logistic Regression Fit") +
  theme_minimal()

model <- glm(cbind(ALT, REF) ~ AI_August, data = df, family = binomial)

summary(model)
# Add predicted values to data frame
df$predicted <- predict(model, type = "response")

ggplot(df, aes(x = AI_August, y = freq)) +
  geom_point(size = 2) +
  geom_line(aes(y = predicted), color = "blue", size = 1) +
  labs(x = "AI (August)", y = "Allele Frequency", title = "Logistic Regression Fit") +
  theme_minimal()

## do individual based
can <- read.table(file = "All_final_Cantabrians.txt", sep="\t", header = F)
env<-read.table(file = "spain_present_day_climate_data_of_all_coordinates.csv", sep=",", header = T)
geno <- read.table(file = "geno_nac055_CAN.txt", sep="\t", header = F)
geno <- read.table(file = "geno_nac055_10078172_CAN.txt", sep="\t", header = F)
mds <- read.table(file = "plink_tab.mds", sep="\t", header = T)



env <- env[env$VCF.BGI.code %in% geno$V1,]

env<- env[order(env$VCF.BGI.code),]
mds <- mds[order(mds$FID),]
geno <- geno[order(geno$V1),]
can <- can[order(can$V1),]

nrow(geno)
nrow(env)
nrow(mds)
nrow(can)

df <- cbind(geno,can$V6, env$ai_Jul/10000, mds$C1,mds$C2,mds$C3,mds$C4)
colnames(df) <- c("geno", "dosage", "pop", "AI_july", "C1", "C2", "C3", "C4")

df$group <- "W-CAN1"
df$group[df$pop %in% c("ES04", "ES05", "ES23")] <- "W-CAN2"
df$group[df$pop %in% c("ES17", "ES08", "ES09", "ES24","ES25")] <- "CE-CAN"

model<-glm(cbind(dosage, 2 - dosage) ~ AI_july + C1 + C2 + C3 + C4,
    family = binomial,
    data = df)
m1_all<-summary(model)


# Generate new data for prediction
newdata <- data.frame(
  AI_july = seq(min(df$AI_july), max(df$AI_july), length.out = 100),
  C1 = 0, C2 = 0, C3 = 0, C4 = 0  # hold structure constant
)

# Predict fitted probabilities (ALT allele freq)
newdata$fit <- predict(model, newdata, type = "response")

# Plot


ggplot(df, aes(x = AI_july, y = dosage/2, color = group, shape = group)) +
  geom_point(alpha = 0.5, size = 4) +
  geom_line(data = newdata, aes(x = AI_july, y = fit), inherit.aes = FALSE, color = "red", size = 1.2) +
  scale_y_continuous(breaks = c(0, 0.5, 1), labels = c(0, 1, 2)) +
  scale_x_continuous(breaks = seq(0.25, 0.4, by = 0.05)) +
  labs(
    x = "AI july",
    y = "Derived allele dosage",
    color="",
    shape=""
  ) +
  geom_segment(y = -0.05, yend = -0.05, x = 0.25, xend = 0.4, color="black") +
  geom_segment(y = 0, yend = 1, x = 0.24175, xend = 0.24175, color="black") +
  scale_color_manual(values = c("#1f78b4", "#e66101", "#33a02c")) +
  scale_shape_manual(values = c(16, 17, 15)) +
  theme(
    axis.text.x = element_text(size = 40),
    axis.text.y = element_text(size = 40),
    panel.background = element_rect(fill = NA, color = "white"),
    axis.title.y = element_text(size = 40, vjust = 3),
    axis.title.x = element_text(size = 40, vjust = -2),
    legend.title = element_text(size = 40),
    legend.key.size = unit(2, "cm"),
    legend.key = element_rect(fill = "white"),
    legend.text = element_text(size = 40),
    axis.ticks.length.y = unit(0.5, "cm"),
    axis.ticks.length.x = unit(0.5, "cm"),
    plot.margin = unit(c(1, 1, 1, 1), "cm"),
    plot.title = element_text(size = 40, hjust = 0.5, vjust = 1.5),
    legend.position = "top"
  )+
  guides(
    color = guide_legend(override.aes = list(size = 10, alpha = 1)),
    shape = guide_legend(override.aes = list(size = 10, alpha = 1))
  )


AI_july_p1 <- ggplot(df, aes(x = AI_july, y = dosage/2, color = group, shape = group)) +
  geom_point(alpha = 0.5, size = 6) +
  geom_line(data = newdata, aes(x = AI_july, y = fit), inherit.aes = FALSE, color = "red", size = 1.2) +
  scale_y_continuous(breaks = c(0, 0.5, 1), labels = c(0, 1, 2)) +
  scale_x_continuous(breaks = seq(0.25, 0.4, by = 0.05)) +
  labs(
    x = "AI july",
    y = "Derived allele dosage",
    color="",
    shape=""
  ) +
  geom_segment(y = -0.05, yend = -0.05, x = 0.25, xend = 0.4, color="black") +
  geom_segment(y = 0, yend = 1, x = 0.24175, xend = 0.24175, color="black") +
  scale_color_manual(values = c("#1f78b4", "#e66101", "#33a02c")) +
  scale_shape_manual(values = c(16, 17, 15)) +
  theme(
    axis.text.x = element_text(size = 40),
    axis.text.y = element_text(size = 40),
    panel.background = element_rect(fill = NA, color = "white"),
    axis.title.y = element_text(size = 40, vjust = 3),
    axis.title.x = element_text(size = 40, vjust = -2),
    legend.title = element_text(size = 40),
    legend.key.size = unit(2, "cm"),
    legend.key = element_rect(fill = "white"),
    legend.text = element_text(size = 40),
    axis.ticks.length.y = unit(0.5, "cm"),
    axis.ticks.length.x = unit(0.5, "cm"),
    plot.margin = unit(c(1, 1, 1, 1), "cm"),
    plot.title = element_text(size = 40, hjust = 0.5, vjust = 1.5),
    legend.position = "none"
  )+
  guides(
    color = guide_legend(override.aes = list(size = 10, alpha = 1)),
    shape = guide_legend(override.aes = list(size = 10, alpha = 1))
  )

## exclude ES10, ES02, ES15
exclude<-can$V1[can$V6 %in% c("ES10", "ES02", "ES15")]

df<- df[!df$geno%in%exclude,]

model<-glm(cbind(dosage, 2 - dosage) ~ AI_july + C1 + C2 + C3 + C4,
           family = binomial,
           data = df)
summary(model)


# Generate new data for prediction
newdata <- data.frame(
  AI_july = seq(min(df$AI_july), max(df$AI_july), length.out = 100),
  C1 = 0, C2 = 0, C3 = 0, C4 = 0  # hold structure constant
)

# Predict fitted probabilities (ALT allele freq)
newdata$fit <- predict(model, newdata, type = "response")

# Plot
AI_july_p2 <- ggplot(df, aes(x = AI_july, y = dosage/2, color = group, shape = group)) +
  geom_point(alpha = 0.5, size = 6) +
  geom_line(data = newdata, aes(x = AI_july, y = fit), inherit.aes = FALSE, color = "red", size = 1.2) +
  scale_y_continuous(breaks = c(0, 0.5, 1), labels = c(0, 1, 2)) +
  scale_x_continuous(breaks = seq(0.25, 0.4, by = 0.05)) +
  labs(
    x = "AI july",
    y = "Derived allele dosage",
    color="",
    shape=""
  ) +
  geom_segment(y = -0.05, yend = -0.05, x = 0.25, xend = 0.4, color="black") +
  geom_segment(y = 0, yend = 1, x = 0.24175, xend = 0.24175, color="black") +
  scale_color_manual(values = c("#1f78b4", "#e66101", "#33a02c")) +
  scale_shape_manual(values = c(16, 17, 15)) +
  theme(
    axis.text.x = element_text(size = 40),
    axis.text.y = element_text(size = 40),
    panel.background = element_rect(fill = NA, color = "white"),
    axis.title.y = element_text(size = 40, vjust = 3),
    axis.title.x = element_text(size = 40, vjust = -2),
    legend.title = element_text(size = 40),
    legend.key.size = unit(2, "cm"),
    legend.key = element_rect(fill = "white"),
    legend.text = element_text(size = 40),
    axis.ticks.length.y = unit(0.5, "cm"),
    axis.ticks.length.x = unit(0.5, "cm"),
    plot.margin = unit(c(1, 1, 1, 1), "cm"),
    plot.title = element_text(size = 40, hjust = 0.5, vjust = 1.5),
    legend.position = "none"
  )+
  guides(
    color = guide_legend(override.aes = list(size = 10, alpha = 1)),
    shape = guide_legend(override.aes = list(size = 10, alpha = 1))
  )

nrow(df)
library(cowplot)
plot_grid(plotlist = list(AI_july_p1,AI_july_p2),align = "vh", nrow = 1)


######
df$C1
model <- lm(AI_july ~ dosage + C1 + C2 + C3 + C4, data = df)
summary(model)

