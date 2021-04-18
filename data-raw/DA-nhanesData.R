#' @importFrom magrittr %>%
#' @importFrom dplyr filter bind_rows arrange mutate
#' @importFrom stats na.omit
#' @importFrom RNHANES nhanes_data_files nhanes_search nhanes_load_data

##Use this block just once to load and prepare the raw nhanes dataset
#load NHANES datasets for body measurements with demographics data
dat_tot <- nhanes_data_files()  #load NHANES data files
dat_BMX <- nhanes_search(dat_tot, "BMX")  #get datasets with BMX (body measurements)

#make new dataframe with required file names and cycles
df <- data.frame(file_name=dat_BMX$data_file_name, cycles=dat_BMX$cycle, stringsAsFactors=FALSE)
df <- df %>% filter(cycles != "2015-2016")  #get rid of problematic "2015-2016" record

#extract desired data
dat_list <- nhanes_load_data(df$file_name, year=df$cycles, demographics = TRUE) #get data list
dat <- bind_rows(dat_list)  #bind to a single dataframe

#extract the desired variables; for reference on variable names go to: https://www.cdc.gov/nchs/nhanes/index.htm
#BMXWT = weight (kg)
#BMXHT = height (cm)
#BMXBMI = body mass index (kg/m^2)
#RIAGENDR = gender (1=male, 2=female)
#RIDAGEYR = age (years)
#RIDAGEMN = age (months)
#RIDRETH2 = race (1=non-Hispanic white, 2=non-Hispanic black, 3=Mexican American, 4=other race, 5=other Hispanic)

dat2 <- data.frame(SEX=dat$RIAGENDR, AGE_YR=dat$RIDAGEYR, AGE_MN=dat$RIDAGEMN, RACE=dat$RIDRETH2, BW=dat$BMXWT, HT=dat$BMXHT, BMI=dat$BMXBMI)
dat2 <- dat2 %>% na.omit()  #get rid of NAs; note: this will get rid of all subjects < 2 yrs of age (no height)
nhanesDataSet <- dat2 %>% mutate(BSA = 0.024265 * BW^0.5378 * HT^0.3964) %>% arrange(AGE_YR) #calculating BSA using Haycock formula
