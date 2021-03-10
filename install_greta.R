install.packages("reticulate",dependencies = TRUE)
install.packages("greta",dependencies = TRUE)
install.packages("causact",dependencies = TRUE)

condaPath = try(reticulate::conda_binary()) #error is OK
if (is.null(condaPath) | class(condaPath) == "try-error") {
  ## install miniconda in default location if possible
  condaInstall = try(reticulate::install_miniconda())
  condaPath = try(reticulate::miniconda_path())
  ## if install still fails install miniconda in C:\miniconda
  if(class(condaInstall) == "try-error") {
    condaPath = file.path("/", "miniconda")
    reticulate::install_miniconda(path = condaPath,
                                  force = TRUE)}
  reticulate::miniconda_update(path = condaPath) #update new install
  newCondaInstall = TRUE  ## flag to update env var if needed
} else {
  reticulate::miniconda_update() ## update existing conda
  newCondaInstall = FALSE
  ### this environment is added automatically by install_miniconda()
  reticulate::conda_create(envname = "r-reticulate") ## add environment
}

#### STEP 3: RESTART R WITH CONDA ENVIRONMENT HINT
#### GIVEN IN EDITED .Rprofile FILE
## capture environment variable for inclusion in .Rprofile
newProfileLines = c() ## initialize config lines for .Rprofile
if(newCondaInstall) {  ## help R find new install location
  newProfileLines = c(
    newProfileLines,
    paste0("Sys.setenv(RETICULATE_MINICONDA_PATH = '",
           condaPath, "')")
  )
  rProfilePath = file.path("~", ".Rprofile")
  profileLines = c()  ## init blank lines
  if (file.exists(rProfilePath)) {
    ## is there a .RProfile file?
    profileLines = readLines(rProfilePath)# get rProfile
  }
  ## add new line to top of file
  profileLines = c(newProfileLines, profileLines)
  writeLines(profileLines, rProfilePath)
  .rs.restartR()  ## run & wait for R session to resume 
  ## (at least 15 seconds) before continuing
}

#### STEP 4:  Update "r-reticulate" CONDA ENVIRONMENT
####          FOR TENSORFLOW
## Install the specific versions of modules
## for the TensorFlow installation via CONDA.
## this line takes a few minutes to execute
reticulate::conda_install(envname = "r-reticulate",
                          packages =
                            c(
                              "python=3.7",
                              "tensorflow=1.14",
                              "pyyaml",
                              "requests",
                              "Pillow",
                              "pip",
                              "numpy=1.16",
                              "h5py=2.8",
                              "tensorflow-probability=0.7"
                            ))

#### STEP 4:  TEST THE INSTALLATION
library(greta)  ## should work now
library(causact)
library(dplyr)
graph = dag_create() %>%
  dag_node("Normal RV",
           rhs =normal(0,10))
graph %>% dag_render()  ## see oval
drawsDF = graph %>% dag_greta() ## see "running X chains..."
drawsDF %>% dagp_plot(densityPlot = TRUE)  ## see plot