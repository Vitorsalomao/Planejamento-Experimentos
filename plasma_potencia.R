
library(dplyr)
library(ggplot2)


rm(list = ls())
plasmac <- c(575, 542, 530, 539, 570, 565, 593,
             590, 579, 610, 600, 651, 610, 637,
             629, 725, 700, 715, 685, 710)
length(plasmac)
potencia <- c(rep(160,5), rep(180,5),
              rep(200,5), rep(220,5))
length(potencia)
pdata <- data.frame(plasmac = plasmac, potencia = as.factor(potencia))
View(pdata)
result <- pdata %>%
  group_by(potencia) %>%
  summarise(mean_plasmac = mean(plasmac))                      

ggplot(result, aes(x = potencia, y = mean_plasmac)) +
  geom_bar(stat = "identity") +
  labs(title = "Média de plasma por Potência",
       x = "Potência",
       y = "Média de Plasma")

gp.means <- with(pdata, tapply(plasmac, potencia, mean))
with(pdata, stripchart(plasmac ~ potencia, vert = T, method = 'overplot',
                       pch = 1, ylab = ''))
stripchart(plasmac ~ potencia, data = pdata, vertical = TRUE, method = 'overplot',
           pch = 1, xlab = 'Potência', ylab = 'Plasma')

##MODELO DE EFEITOS
#Yij = M + T + eji
##MODELO DE MÉDIAS
#Yij = Mi + eij

plasmac.aov <- aov(plasmac~potencia, data = pdata)
summary(plasmac.aov)

##Estimador de T = Yi.(BARRA) - Y..(BARRA)
#Média das classes - Média geral
View(result)
T1 =  result$mean_plasmac[1] - mean(plasmac); T1
T2 =  result$mean_plasmac[2] - mean(plasmac); T2
T3 =  result$mean_plasmac[3] - mean(plasmac); T3
T4 =  result$mean_plasmac[4] - mean(plasmac); T4
#Ou de forma mais simples
T_data = gp.means - mean(plasmac); T_data
#Ou ainda
model.tables(plasmac.aov)
## INTERVALO DE CONFIANÇA
#IC(Mi;1-alpha) = ( Yi.(BARRA) -+ t-student(1-alpha/2);N-a*sqrt(QMe/n))

quantil <- qt(.975, 16) #Quantil para 97.5% e 16 GL
qme <- sigma(plasmac.aov)^2;qme # QUADRADO MÉDIO DO RESÍDUO
li =  gp.means - quantil*sqrt(qme/5);li #limite inferior
ls =  gp.means + quantil*sqrt(qme/5);ls #limite superior
confint(plasmac.aov) # ATENÇÃO: Intervalo de confiança da Média i com a Média 1

##SUPOSIÇÕES DO MODELO
#Avaliando a Normalidade com o Q-Q Plot com os resíduos
qqnorm(plasmac.aov$residuals)
qqline(plasmac.aov$residuals)
#Teste Shapiro Wilk
shapiro.test(plasmac)
#Avaliando a Homocedasticidade
#Teste de BARTLETT 
bartlett.test(plasmac ~ potencia, data = pdata)
#Teste de Levene
library(car)
leveneTest(plasmac.aov)
#Teste Durbin–Watson
set.seed(474028)
durbinWatsonTest(plasmac.aov) # Ele utiliza Bootstrap, então é necessário fixar uma semente

# Avaliando a comparação entre as médias
##H0: mi=mj
##H1: mi != mj
pairwise.t.test(pdata$plasmac, pdata$potencia, p.adjust.method = 'bonf')

#Intervalo e coniança para a diferença das médias
#Teste Tukey 
Tk <- TukeyHSD(plasmac.aov);Tk 
plot(Tk, las = 1) #Gráfico para todas diferenças das médias; Se o intevalo de confiança não contém o zero rejeita-se H0

#ANALISANDO O CONTRASTE
stc <- ScheffeTest(plasmac.aov);stc
stc <- ScheffeTest(plasmac.aov, contrasts = matrix(c(1,-1,0,0,1,-0.5, -0.5,0, 0, 0, 1, -1), ncol=3));stc
