listlibraries = c("TSstudio","ggplot2","dplyr","xts","forecast","tseries")

loadpkg = function(x){
    options(warn=-1)
    for(i in x){
        tryCatch({
            library(i,character.only = TRUE)
        }, error = function(e){
            cat("Installing ",i)
            install.packages(i,dependencies = TRUE,repos='http://cran.us.r-project.org')
            library(i,character.only = TRUE)
        }
        )}
    options(warn=0)
    cat("\014") 
    cat("Successfully install and load all packages along with its dependencies")
}

loadpkg(listlibraries)

test_size <<- 0.15
glob.alpha <<- 0.05
set.seed(28062019)

df = read.csv("datameetupr (1).csv",header=T,sep=",")
colnames(df) =  c("date","active_users")

head(df)

df = xts(df$active_users, order.by=as.Date(df$date,"%Y-%m-%d"))
head(df)

#Read All
df_all = read.csv("datameetupr (2).csv",header=T,sep=",")
head(df_all)

colnames(df_all) = c("ds","y")
df_all$ds = as.POSIXct((df_all$ds),tz = "GMT",format="%Y-%m-%d")
head(df_all)

#Data Training
split_boundary = (1-test_size)*dim(df_all)[1]
df_train = df_all[1:split_boundary,]
tail(df_train)
print(dim(df_train)[1])

#Data Testing
df_test = df_all[(split_boundary+1):dim(df_all)[1],]
head(df_test)
print(dim(df_test)[1])

fit_arima = auto.arima(df_train$y)

prediction_arima = forecast(fit_arima, h=134)
plot(prediction_arima)

accuracy(fit_arima)
