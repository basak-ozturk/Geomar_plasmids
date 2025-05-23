pacman::p_load(tidyverse)

testset1 <- c("Meier","Mayer","Maier","Meyer","Mayr","Maya","Mayor")
# find all variations of the name "Meier"

testset2 <- c("weight_mm","height_cm","age_yr","temp_c")
#replace _ with space
#replace _ with space and add unit in brackets


testset3 <- c("1980_12_30","13.04.2005", "2005/04/25","24121990")
# transform into YYYY-MM-DD

testset4 <- c("pw2000","That1sb3tt3r","M@kesSense?","NoDigits@this1")
# test pwd strength, rules: Upper, lower, special char, number, min 8 char long


str_detect(string = testset1, pattern = "M[ea][iy]e?r$")

str_replace(string=testset2, 
            pattern = "_", 
            replacement = " ")

new_testset2<-str_replace(string=testset2, 
            pattern = "_", 
            replacement = " \\[")

print(new_testset2)

final_testset2 <- str_replace(string = new_testset2, 
                              pattern = "$", 
                              replacement = "]")

print(final_testset2)

