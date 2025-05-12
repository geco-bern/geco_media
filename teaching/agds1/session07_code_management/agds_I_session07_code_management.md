# PLAYBOOK 2025-03-31 --- Presentation

## In Action: live presentation of git workflow

We will:

1. start a new GitHub repo
2. create a connected RStudio project
3. write the code for our data analysis 
4. New R script in /analysis/01_plot.R and write code and 
5. Push it to github

### How to make a new GitHub repo:

>  github.com, new (without template)
> agds_example_git

### How to make a RStudio git project:

RStudio: New Project => Version Control => git => `git@github.com:fabern/agds_example_git.git`

### Do some data science work

#### 1) New data set in /data/ 

```
dir.create("data")
httr::GET("https://raw.githubusercontent.com/geco-bern/agds_book/refs/heads/main/book/data/FLX_CH-Lae_FLUXNET2015_FULLSET_HH_2004-2006_CLEAN.csv",
	httr::write_disk(path = "data/FLX_CH-Lae_FLUXNET2015_FULLSET_HH_2004-2006_CLEAN.csv",overwrite = TRUE))
```

##### 2) New R script in `/analysis/01_plot.R`

```R
library(readr)
library(ggplot2)

half_hourly_fluxes <- readr::read_csv("../data/FLX_CH-Lae_FLUXNET2015_FULLSET_HH_2004-2006_CLEAN.csv")

# prepare plot data
plot_data <- half_hourly_fluxes |> 
  dplyr::slice(24000:25000)

# plot figure
plotme <- ggplot(
    data = plot_data,
    aes(x = TIMESTAMP_START, y = GPP_NT_VUT_REF)) +
  geom_line() +
  labs(title = "Gross primary productivity", 
       subtitle = "Site: CH-Lae",
       x = "Time", 
       y = expression(paste("GPP (", mu,"mol CO"[2], " m"^-2, "s"^-1, ")"))) +
  theme_classic()
plotme
```
> Commit with RStudio diff: "Generate plot of GPP May and June 2005"
> Show History
> Show that Git panel is now emtpy again
> Then go on and do some more work

##### 2b) Generate output and more statistics:

```R
ggsave(plot = plotme, "analysis/figures/01_plot.png")

# output statistics
library(lubridate)
library(dplyr)
stats <- half_hourly_fluxes |> 
  mutate(year = lubridate::year(TIMESTAMP_START),
         month = lubridate::month(TIMESTAMP_START)) |>
  group_by(year, month) |> summarise(N = n(),
                                     mean_GPP = mean(GPP_NT_VUT_REF))
stats

dir.create("analysis/statistics")
readr::write_csv(stats, "analysis/statistics/01_plot.csv")
```
> Run code. (generates png and csv file)
> Commit only Rscript changes: "Compute monthly statistics and save figures and statistics"
> For the moment don't consider output files
> Continue working

##### 2c) Re-organize code: move `library(lubrdiate)`, `library(dplyr)` to top of file

> commit: "Reorganize code"

> amend .gitignore

```
.Rproj.user
.Rhistory
.Rdata
.httr-oauth
.DS_Store
.quarto

# Output files from analysis
analysis/figures/*
analysis/statistics/*
```

> commit: "Add outputs to .gitingore"

So far everything has been local on my computer.

#### 3) Create github project and push it to github

```R
usethis::use_github()
# note that this is different from usethis::use_git()

# might need to do:
usethis::create_github_token()
# "usethis RStudio T16" and then copy the created token
gitcreds::gitcreds_set()
```

> go to https://github.com/fabern/agds_example_git_started
> show that files are now here
> show that outputs were not pushed
> show that collaborator can clone this repository

### **Collaborate:** 

- anybody: by copying 'forking' your repository, and then suggesting changes to the original copy ('pull request')
- a collaborator: directly in your repository if you allowed them as collaborator (use branches.)

Once you start working on the same project, you will need to start
`pull`ing the changes from GitHub before you start working
and 
`push`ing them after you have finished.

**Conflicts:** git is line based, so if you and your collaborator don't modify the same lines in a file at the same time, it will simply work.
However, if both of you modified the same files so-called `merge conflicts` can happen. These need to be resolved manually.

### Report exercise:

- can be based on the exercise (do that first)
- collaboration with a partner: make groups of 2
- do a fork, a pull request
