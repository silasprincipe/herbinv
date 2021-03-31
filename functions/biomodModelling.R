###### Reef-builders modelling #######
#Silas C. Principe - 2020
#silasprincipe@yahoo.com.br

# Function to execute the modeling using BIOMOD2 ----

biomodModelling <- function(species, round.code){
        
        #Load libraries ----
        library(raster)
        library(tidyverse)
        library(biomod2)
        
        # Set seed for replicability
        set.seed(2932)
        
        # Load environmental layers (present and future) ----
        source("functions/varload.r")
        env <- var.load(folder = "crop_layers/", layers = "env_layers.txt")
        #env <- dropLayer(env, "BO21_lightbotmax_bdmax")
        
        r26 <-
                var.load(folder = "proj_layers/ssp126/",
                        layers = "env_layers.txt")

        r45 <-
                var.load(folder = "proj_layers/ssp245/",
                        layers = "env_layers.txt")

        r85 <-
                var.load(folder = "proj_layers/ssp585/",
                        layers = "env_layers.txt")
        
        
        # Define modelling parameters ----
        
        # Algorithms to use
        models <- c("RF", "GBM", "GLM")
        
        # Variable importance permutation number
        vimp <- 20
        
        # TSS threshold
        tss.th <- 0.7
        
        # Run code (will be the folder name)
        pn <- paste(species, round.code, sep = "_")
        
        # Set model options (we use the standard of BIOMOD2)
        biomod.options <- BIOMOD_ModelingOptions()
        
        # Register the start time of modelling
        ptm.t <- proc.time()
        
        print(paste(species, "model started at", Sys.time(),
                    "! Run code", pn, sep = " "))
        
        # Load species data ----
        pts <-
                read.csv(paste("data/", species, "/",
                               species, "_pts.csv", sep = ""))
        set <-
                read.csv(paste("data/", species, "/",
                               species, "_set.csv", sep = ""))
        dsp <-
                read.csv(paste("data/", species, "/",
                               species, "_dsp.csv", sep = ""))
        
        
        # Load split table (block cross validation)
        split.table <-
                read.csv(paste0("data/", species, "/",
                                species, "_splittable.csv"))
        split.table$RUN6 <- TRUE
        split.table <- as.matrix(split.table)
        
        # Create the data object to BIOMOD ----
        biomod.data <- BIOMOD_FormatingData(
                resp.var = dsp,
                expl.var = env,
                resp.xy = pts,
                resp.name = species,
                PA.strategy = 'user.defined',
                PA.table = set
        )
        
        
        # Run modeling ----
        biomod.model.out <- BIOMOD_Modeling(
                biomod.data,
                models = models,
                models.options = biomod.options,
                Prevalence = 0.5,
                VarImport = vimp,
                models.eval.meth = c('TSS', 'ROC', 'KAPPA'),
                SaveObj = TRUE,
                DataSplitTable = split.table,
                modeling.id = pn
        )
        
        # Get models evaluation
        biomod.model.eval <-
                get_evaluations(biomod.model.out, as.data.frame = T)
        
        # Select the models with TSS score higher than the threshold
        sel.models <-
                separate(biomod.model.eval, Model.name, c("model", "run", "pa"),
                         sep = "_")
        sel.models <-
                sel.models %>% filter(Eval.metric == "TSS") %>%
                filter(run != "RUN6") %>% group_by(model) %>%
                summarise(mean = mean(Testing.data))
        sel.models <- sel.models %>% filter(mean >= tss.th)
        
        # Save selected models
        if (dir.exists(paste(species, "/evals/", sep = "")) == F) {
                dir.create(paste(species, "/evals/", sep = ""), recursive = T)
        }
        
        write.csv(sel.models,
                  file = paste(species, "/evals/", pn,
                               "_selmodelstotal.csv", sep = ""))
        
        # Get the names of selected models for projection
        sel.models <- paste(species, "_PA1_RUN6_", sel.models$model, sep = "")
        
        
        # Projection in current conditions ----
        biomod.proj <- BIOMOD_Projection(
                modeling.output = biomod.model.out,
                new.env = env,
                proj.name = paste("current", pn, sep = "_"),
                selected.models = sel.models,
                binary.meth = 'TSS',
                compress = F,
                build.clamping.mask = F,
                do.stack = F,
                output.format = '.img',
                on_0_1000 = T
        )
        
        # Ensemble modelling and projection in current conditions ----
        # Ensemble of models
        biomod.ensemble <- BIOMOD_EnsembleModeling(
                modeling.output = biomod.model.out,
                chosen.models = sel.models,
                em.by = 'all',
                eval.metric = c('TSS'),
                eval.metric.quality.threshold = tss.th,
                prob.mean = F,
                prob.cv = F,
                prob.ci = F,
                prob.ci.alpha = 0.05,
                prob.median = F,
                committee.averaging = T,
                prob.mean.weight = T,
                prob.mean.weight.decay = 'proportional',
                VarImport = vimp
        )
        
        
        # Ensemble Forecasting
        biomod.ens.proj <-
                BIOMOD_EnsembleForecasting(
                        projection.output = biomod.proj,
                        EM.output = biomod.ensemble,
                        output.format = '.img',
                        on_0_1000 = T
                )
        
        # Projection in future conditions ----
        # Create a function that make ensemble and projection in the future
        futureProj <- function(scenario){

                if(scenario == "ssp126"){
                        fut.env <- r26
                }
                if(scenario == "ssp245"){
                        fut.env <- r45
                }
                if(scenario == "ssp585"){
                        fut.env <- r85
                }


                # Ensemble projection in future condition
                biomod.proj.fut <- BIOMOD_Projection(
                        modeling.output = biomod.model.out,
                        new.env = fut.env,
                        proj.name = paste(scenario, pn, sep = "_"),
                        selected.models = sel.models,
                        binary.meth = 'TSS',
                        compress = F,
                        build.clamping.mask = F,
                        do.stack = F,
                        output.format = '.img',
                        on_0_1000 = T
                )

                # Ensemble Forecasting
                biomod.ens.proj <-
                        BIOMOD_EnsembleForecasting(
                                projection.output = biomod.proj.fut,
                                EM.output = biomod.ensemble,
                                output.format = '.img',
                                on_0_1000 = T
                        )
        }

        # # Apply for all scenarios
        lapply(c("ssp126", "ssp245", "ssp585"), futureProj)
        
        # Save model evaluations in CSV format ----
        
        # Model evaluations
        biomod.model.eval <-
                get_evaluations(biomod.model.out, as.data.frame = T)
        write.csv(
                biomod.model.eval,
                paste(species, "/", "evals/model_eval_", pn, ".csv", sep = ""),
                row.names = F
        )
        
        # Variable importances
        biomod.var.imp <-
                as.data.frame(get_variables_importance(biomod.model.out,
                                                       as.data.frame = T))
        write.csv(
                biomod.var.imp,
                paste(species, "/", "evals/var_imp_", pn, ".csv", sep = ""),
                row.names = T
        )
        
        # Ensemble variable importances
        biomod.var.imp <-
                get_variables_importance(biomod.ensemble)
        
        write.csv(
                biomod.var.imp[,,1],
                paste(species, "/", "evals/ensemble_wm_var_imp_", pn,
                      ".csv", sep = ""),
                row.names = T
        )
        
        # Ensemble evaluations
        ensemble.evaluations <-
                get_evaluations(biomod.ensemble, as.data.frame = T)
        write.csv(
                ensemble.evaluations,
                paste(species, "/", "evals/ensemble_eval_", pn,
                      ".csv", sep = ""),
                row.names = F
        )
        
        ### End of modelling code
        
        timed.t <- proc.time() - ptm.t
        t.min.t <- floor(timed.t[3] / 60)
        t.sec.t <- timed.t[3] - (t.min.t * 60)
        
        cat(paste("\n", "Model", pn, "done!", "\n", sep = " "))
        
        cat(species,
            "modelling finished in",
            t.min.t,
            "minutes",
            round(t.sec.t, 1),
            "seconds.",
            "Round code:",
            round.code)
        
        gc()
        
}

###END