################################################################################
# Aula sobre a biblioteca spatstat e splancs:
#################################################################################

#Carregando o pacote
library(spatstat)
library(splancs)                                         


# O spatstat e splancs são pacotes para trabalhar com dados espaciais ligados a 
# processos pontuais.


################################################################################
# Visualização de processos pontuais:
################################################################################

# Para a biblioteca spatstat #

# Uma demonstração do pacote pode ser vista, da seguinte maneira :
# demo(spatstat)


# A biblioteca trabalha com objetos da classe "ppp". Se tivermos dados 
# e quisermos utilizar esta função, temos que transformá-los em objeto "ppp".

# Exemplo: criando coordenadas fictícias
x <- runif(100)
y <- runif(100)
exemplo <- ppp(x, y)
str(exemplo)
plot(exemplo,pch=23,bg=2)
quadrantes1 <- quadratcount(exemplo,nx=4,ny=4)    
quadrantes1 

par(mar=c(0,0,0,0))
plot(exemplo,pch="*",main="",fg="red")
plot(quadrantes1,add=TRUE,col="red", cex=1.5, lty=2,border="red")


# O default é criar em um quadrado, mas pode ser criado em figuras poligonais
# diversas
exemplo2 <- ppp(x, y, poly=list(x=c(0,2,0), y=c(0,0,2)))
plot(exemplo2)
quadrantes2 <- quadratcount(exemplo2,nx=4,ny=4)    
quadrantes2 

par(mar=c(0,0,0,0))
plot(exemplo2,pch="*",main="",fg="red")
plot(quadrantes2,add=TRUE,col="red", cex=1.5, lty=2,border="red")

# Outro comando para criar processos pontuais uniformes
pontos <- runifpoint(100,win=owin(c(0,1),c(0,1)))
plot(pontos,pch=23,bg=2)
 
# Já existem alguns arquivos pré-criados.
# Entre eles:
# amacrine
# ants
# nztrees
# redwood

# Para carregá-los basta fazer:
data(cells) # localização do centroide de células vista em um microscópio #
str(cells)
plot(cells)

data(redwood) # Localização de árvores de madeira vermelha na California #
redwood
plot(redwood)


data(amacrine) # um padrão de ponto de células amácrinas deslocadas na retina de 
               # um coelho. 152 células "ativas" e 142 células "desativadas" 
               # em uma estrutura de amostragem retangular. #
str(amacrine)
plot(amacrine)

data(ants) # ninhos de duas espécies de formigas em um local no norte da Grécia #
str(ants)
plot(ants)


# Para a biblioteca splancs #

# A biblioteca trabalha com objetos da classe "pts". Se tivermos dados 
# e quisermos utilizar esta função, temos que transformá-los em objeto "pts".

#Atenção mudar o diretório para onde foram criados os arquivos !!!!

coord <- read.table("coord.txt")
polig <- read.table("pinn.txt")

exemplo3 <- splancs::as.points(coord[,1],coord[,2])

plot(polig, asp=1, type="n")
splancs::polymap(polig, add=TRUE, col=3)
splancs::pointmap(exemplo3,add=TRUE, bg=2, pch = 21)


#Também existem alguns dados já existentes no pacote

data(uganda)   #localização de crateras vulcânicas em Uganda#

par(mar=c(4,4,0.5,0.5))
plot(uganda$poly, asp=1, type="n")
points(as.points(uganda),pch=16)
#points(as.points(csr(uganda$poly,120)),pch="*",lwd=10,col=2) # gerando um processo CSR
polymap(uganda$poly,add=TRUE)


################################################################################
# Explorando e modelando dados de processos pontuais:
################################################################################


#########################
# Método dos quadrantes #
#########################

# Na biblioteca spatstat
data(bei)
plot(bei)
quadrantes1 <- quadratcount(bei,nx=5,ny=4)    
quadrantes1 


plot(bei,pch="+")
plot(quadrantes1,add=TRUE, col="red", cex=1.5, lty=2)

t_quadrantes <- quadrat.test(bei,nx=5,ny=4)    
t_quadrantes
qchisq(0.95,19)



quadrantes2 <- quadratcount(cells,nx=5,ny=3) 

plot(cells,pch="+")
plot(quadrantes2,add=TRUE, col="red", cex=1.5, lty=2)

t_quadrantes <- quadrat.test(cells,5,3) 
t_quadrantes
qchisq(0.95,14)

# Índices
media1 <- mean(quadrantes1)
variancia1 <- sum((quadrantes1 - media1)^2/19)

media2 <- mean(quadrantes2)
variancia2 <- sum((quadrantes2 - media2)^2/14)

# ICS:
(variancia1/media1) - 1

# ICS:
(variancia2/media2) - 1               



# Intensidade estimada lambda_i
intensity(quadrantes1)
par(mar=c(0,0,0,2))
plot(intensity(quadrantes1, image=TRUE), main="", 
     col=terrain.colors(256))



#######################################
# Estimador de intensidade via KERNEL #
#######################################

# Para o pacote spatstat


data(redwood) # Localização de árvores de madeira vermelha na California #
plot(redwood,pch="+")

#Estimador de Kernel: Exemplo 1
intensidade <- density.ppp(redwood,sigma=0.05)  #sigma é o raio e podemos colocar efeito borda (edge=TRUE)
plot(intensidade)

intensidade <- density.ppp(redwood,sigma=0.1)  
plot(intensidade)

intensidade <- density.ppp(redwood,sigma=0.5)  
plot(intensidade)

#Exemplo 2:
intensidade2 <- density.ppp(cells,sigma=0.05,edge=TRUE)  
plot(intensidade2)

intensidade2 <- density.ppp(cells,sigma=0.1,edge=TRUE)  
plot(intensidade2)

# Para o pacote splancs

data(bodmin) # Localizações de 35 toras de granito em Bodmin Moor (charneca de granito no nordeste da Cornualha, Inglaterra) 
plot(bodmin$poly, asp=1, type="n")
points(as.points(bodmin),pch=16)
polymap(bodmin$poly,add=TRUE)

image(kernel2d(as.points(bodmin), bodmin$poly, h0=2, nx=200, ny=200),
add=TRUE, col=terrain.colors(20))
pointmap(as.points(bodmin), add=TRUE,pch="*",col=2)
polymap(bodmin$poly, add=TRUE)


data(uganda)
plot(uganda$poly, asp=1, type="n")
image(kernel2d(as.points(uganda), uganda$poly, h0=100, nx=200, ny=200),
add=TRUE, col=heat.colors(20))
pointmap(as.points(uganda), add=TRUE,pch="*")
polymap(uganda$poly, add=TRUE)


plot(uganda$poly, asp=1, type="n")
image(kernel2d(as.points(uganda), uganda$poly, h0=220, nx=200, ny=200),
add=TRUE, col=heat.colors(20))
pointmap(as.points(uganda), add=TRUE,pch="*")
polymap(uganda$poly, add=TRUE)

plot(uganda$poly, asp=1, type="n")
image(kernel2d(as.points(uganda), uganda$poly, h0=500, nx=200, ny=200),
add=TRUE, col=heat.colors(20))
pointmap(as.points(uganda), add=TRUE,pch="*")
polymap(uganda$poly, add=TRUE)


############################
#   Vizinho mais próximo   #
############################

# Para o pacote splancs                            

dists <- nndistG(as.points(uganda))
dists                                                     
maximo <- max(dists$dists)
minimo <- min(dists$dists)
w <- seq(minimo,maximo,by=((maximo - minimo)/100))

# Função G

## o comando Ghat recebe 2 parâmetros, o primeiro o objeto pontos, 
## o segundo um vetor com a lista de distâncias onde se estima a função G

funcao_G <- Ghat(as.points(uganda),w)
plot(w,funcao_G,type="l",xlab="distâncias", ylab="G estimada")

# Função F
dists2 <- nndistF(as.points(csr(uganda$poly, length(uganda$x))), as.points(uganda))
maximo <- max(dists2)
minimo <- min(dists2)
x <- seq(minimo,maximo,by=((maximo - minimo)/100))
funcao_F <- Fhat(as.points(csr(uganda$poly, length(uganda$x))), as.points(uganda),x)
plot(x,funcao_F,type="l", xlab="distâncias", ylab="F estimada")


#Função G x F
plot(funcao_G, funcao_F, type="l", xlab="G estimada", ylab="F estimada")
abline(0,1,lty=2,col=2)



############################
#      Função K            #
############################

# Para a biblioteca splancs:

h <- seq(0,800,10)
funcaoK<- khat(as.points(uganda), uganda$poly, h)
funcaoL<- sqrt(funcaoK/pi) - h

m <- 99
funcaoK.i<-matrix(NA,length(h),m)
funcaoL.i<-matrix(NA,length(h),m)
for (i in 1:m) {
  y<-as.points(csr(uganda$poly, length(uganda$x)))
  funcaoK.i[,i]<- khat(y, uganda$poly, h)
  funcaoL.i[,i]<- sqrt(funcaoK.i[,i]/pi) - h
}
env.sup<-0
for (i in 1:length(h)) {
  env.sup[i]<-max(funcaoL.i[i,])
}
env.inf<-0
for (i in 1:length(h)) {
  env.inf[i]<-min(funcaoL.i[i,])
}

par(mar=c(4.5,4.5,0.5,0.5))
plot(h,funcaoL, type="l", xlab="Distâncias", lwd=2, ylab=expression(hat(L)(h)),ylim=c(-30,60))
lines(h,rep(0,length(h)),lty=2,col=2,lwd=2)
lines(h,env.sup,lty=3,col=4,lwd=2)
lines(h,env.inf,lty=3,col=4,lwd=2)

