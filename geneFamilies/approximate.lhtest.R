library("extRemes")
#The test statistic is calculated as a chi-squared distribution with Degrees of freedom defined as: \[DF={p_1}-{p_2}\] Or, the difference of degrees of freedom from both models.

# s1 should be the model with fewer parameters

#calculate observed likelihood ratio
# Score 4p Lambda: 26377.1
# Score 5p Lambda: 26365.1

s1<-26377.1
s2<-26365.1
lr.test(s1, s2, alpha = 0.05, df = 2)


#calculate observed likelihood ratio
# Score 5p Lambda: 26365.1
# Score 6p Lambda: 26346.4

s1<-26365.1
s2<-26346.4
lr.test(s1, s2, alpha = 0.05, df = 2)


#calculate observed likelihood ratio
# Score 4p Lambda: 26377.1
# Score 5p Lambda: 26365.1
s1<-26377.1
s2<-26365.1
lr.test(s1, s2, alpha = 0.05, df = 2)

#calculate observed likelihood ratio
# Score 6p Lambda: 26346.4
# Score 7p Lambda: 26365.1
s1<-26346.4
s2<-26765.8
lr.test(s1, s2, alpha = 0.05, df = 2)


26365.1


