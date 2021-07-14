# Load libraries ----
library(caret)
library(raster)
library(biomod2)
library(randomForest)
library(dismo)

# Get BIOMOD2 standard parameters ----
std.par <- BIOMOD_ModelingOptions()
std.par@GBM
std.par@RF

# Load data ----
# Establish target species list
sp.list <- c("lyva", "eclu", "trve")

# Load environmental data
source("functions/varload.r")
env <- var.load(folder = "crop_layers/", layers = "env_layers.txt")

# Load species data
sp.data <- lapply(sp.list,
                  function(species){
                          
                          # Load presence points
                          pts <-
                                  read.csv(paste("data/", species, "/",
                                                 species, "_pts.csv", sep = ""))
                          set <-
                                  read.csv(paste("data/", species, "/",
                                                 species, "_set.csv", sep = ""))
                          dsp <-
                                  read.csv(paste("data/", species, "/",
                                                 species, "_dsp.csv", sep = ""))
                          
                          # Convert NA to 0 (we assume pseudo-abs as absences)
                          dsp[is.na(dsp)] <- 0
                          
                          # Extract environmental data for each point
                          full.data <- cbind(species = dsp, raster::extract(env, pts))
                          
                          # Gets the cross-validation index (block cross-val)
                          foldindex <- read.csv("data/lyva/lyva_foldindex.csv")
                          foldindex <- foldindex$x
                          
                          # Create a list to be used in the tunning
                          indexlist <- list(seq(1:3256)[foldindex != 1],
                                            seq(1:3256)[foldindex != 2],
                                            seq(1:3256)[foldindex != 3],
                                            seq(1:3256)[foldindex != 4],
                                            seq(1:3256)[foldindex != 5])
                          
                          # Put everything in a list
                          data <- list(full.data, indexlist)
                          names (data) <- c("data", "index")
                          data
                  })

names(sp.data) <- sp.list

# Train GBM models -----
# Establish control parameters
fit.control <- trainControl(
        method = "cv",
        number = 5,
        index = sp.data[[i]]$index)

# Establish tuning parameters grid
gbm.grid <-  expand.grid(interaction.depth = c(1, 5, 7, 9), 
                        n.trees = seq(1, 20, 0.5)*1000, 
                        shrinkage = seq(0.001, 0.01, 0.005),
                        n.minobsinnode = 5)

# Tune the model
gbm.fit <- train(dsp ~ ., data = sp.data[[i]]$data, 
                 method = "gbm", 
                 trControl = fit.control, 
                 verbose = FALSE, 
                 ## Now specify the exact models 
                 ## to evaluate:
                 tuneGrid = gbm.grid)
gbm.fit


ldata <- sp.data$lyva$data


gbm.tuning <- gbm.step(ldata,
                       gbm.x = 2:8,
                       gbm.y = 1,
                       tree.complexity = 5,
                       learning.rate = 0.002, 
                       n.folds = 4, 
                       max.trees = 20000,
                       bag.fraction = 0.75)

gbm.plot(gbm.tuning)


gbm.prediction <- predict(env, gbm.tuning)

response(gbm.tuning)


#Possibilidade: usar o gbm.step do pacote dismo. Verificar.

# Train randomForest models
# Establish tuning parameters grid
#rf.grid <-  expand.grid(mtry = )

# Tune the model
rf.fit <- tuneRF(ldata[, 2:8], ldata[, 1], ntreeTry = 500)


rf.fit <- randomForest(x = ldata[, 2:8], y = ldata[, 1])


rfTuning <- function(index.list, data, ntree, mtry, nsize){
        
        rmse <- rep(NA, 5)
        
        for (i in 1:5) {
                
                train <- data[index.list[[i]], ]
                test <- data[-index.list[[i]], ]
                
                rf <- randomForest(x = train[, 2:8], y = train[, 1],
                                   ntree = ntree,
                                   mtry = mtry,
                                   nodesize = nsize)
                pred <- predict(rf, newdata = test[, 2:8])
                pred <- round(pred)
                
                RMSE <- function(m, o){
                        sqrt(mean((m - o)^2))
                }
                
                rmse[i] <- RMSE(pred, test[,1])
                
                rm(rf)
        }
        
        final.rmse <- mean(rmse)
        
        final.rmse
}


rf.grid <- expand.grid(ntree = seq(300, 600, 100),
                       mtry = 1:4,
                       nsize = c(1,2,5))

results <- rep(NA, nrow(rf.grid))

pb <- progress_bar$new(total = nrow(rf.grid))

for (j in 1:nrow(rf.grid)) {
        
        results[j] <- rfTuning(sp.data$lyva$index, sp.data$lyva$data,
                               ntree = rf.grid[j,1], mtry = rf.grid[j,2],
                               nsize = rf.grid[j,3])
        pb$tick()
        
        if (j == nrow(rf.grid)) {
                
                bmodel <- rf.grid[which(results == min(results)),]
                
                cat("Best model was with ntree:", bmodel$ntre, "mtry:",
                    bmodel$mtry, "and nodesize:", bmodel$nsize)
        }
        
        
}

#Possibilidade: usar o gbm.step do pacote dismo. Verificar.

rf1 <- randomForest(x = sp.data$lyva$data[, 2:8], y = sp.data$lyva$data[, 1],
                    ntree = 400,
                    mtry = 1,
                    nodesize = 1)

rf2 <- randomForest(x = sp.data$lyva$data[, 2:8], y = sp.data$lyva$data[, 1],
                    ntree = 500,
                    nodesize = 5)


rf3 <- randomForest(x = sp.data$lyva$data[, 2:8], y = sp.data$lyva$data[, 1])



p1 <- predict(env, rf1)
p2 <- predict(env, rf2)
p3 <- predict(env, rf3)


par(mfrow = c(1,3))
plot(p1)
plot(p2)
plot(p3)


# Train GLM models -----
# Establish control parameters
fit.control <- trainControl(
        method = "cv",
        number = 5,
        index = sp.data[[i]]$index)

# Establish tuning parameters grid
gbm.grid <-  expand.grid(interaction.depth = c(1, 5, 7, 9), 
                         n.trees = seq(1, 20, 0.5)*1000, 
                         shrinkage = seq(0.001, 0.01, 0.005),
                         n.minobsinnode = 5)

# Tune the model
glm.fit <- glm(dsp ~ 1 + I(BO21_salinitymin_ss^2) + BO21_tempmean_ss + 
                       I(BO21_tempmean_ss^2) + BO_ph + BO21_lightbotmax_bdmax + 
                       I(BO_ph^2) + I(BO21_silicatemean_ss^2) + BO21_silicatemean_ss + 
                       BO21_chlomean_ss + BO21_dissoxmin_ss,
               data = sp.data[[i]]$data)
gbm.fit


makeFormula("dsp", sp.data[[i]]$data[,2:8], type = "quadratic", interaction.level = 0)








head(ldata)

plot(ldata$BO21_salinitymin_ss, ldata$dsp)

tdata <- ldata %>% mutate("BO21_tempmean_ss" = round(BO21_dissoxmin_ss)) %>%
        group_by(BO21_tempmean_ss) %>%
        summarise("abundance" = sum(dsp))


plot(tdata$BO21_tempmean_ss, tdata$abundance)

ldata$dsp <- as.factor(ldata$dsp)

glm1 <- glm(dsp ~ BO21_tempmean_ss + BO21_salinitymin_ss + BO21_lightbotmax_bdmax +
                    I(BO_ph^2) + BO21_chlomean_ss + BO21_silicatemean_ss + BO21_dissoxmin_ss +
                    BO21_tempmean_ss * BO_ph,
            data = ldata, family = "binomial")
glm1

pchisq(726.1, 3252, lower.tail = F)

BIOMOD_Modeling()

summary(glm1)

head(ldata)

nam <- colnames(ldata)[2:8]

par(mfrow = c(2,2))

for (i in nam) {
        tdata <- ldata %>% mutate("var_new" = round(ldata[,i])) %>%
                group_by(var_new) %>%
                summarise("abundance" = sum(dsp))
        
        colnames(tdata) <- c("var", "abundance")
        
        plot(tdata$var, tdata$abundance)
}

ldata$dsp <- as.factor(ldata$dsp)

glm1 <- glm(dsp ~ BO21_tempmean_ss + I(BO21_tempmean_ss^2) + BO21_salinitymin_ss + BO21_lightbotmax_bdmax +
                    I(BO21_lightbotmax_bdmax^2) +
                    BO_ph + BO21_chlomean_ss + BO21_silicatemean_ss + I(BO21_silicatemean_ss^2) +
                    BO21_dissoxmin_ss + I(BO21_dissoxmin_ss^2) + BO21_tempmean_ss * BO_ph + BO21_tempmean_ss * BO21_salinitymin_ss ,
            data = ldata, family = "binomial")
glm1

library(MASS)

glm.fit <- stepAIC(glm1)

