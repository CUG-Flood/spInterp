#' @title Inverse Distance Weighting interpolation
#' @name IDW
#'
#' @description Apply the Inverse Distance Weighting interpolation
#' method on the irregular points of meteorological observations to a regular grid.
#'
#' @param xy A matrix (N,2) with longitude and latitude of points of data observed 
#' @param z A vector (N) with the values observeds in the points
#' @param xrange A vector with extremes longitudes of grid eg: c(-45,-35)
#' @param yrange A vector with extremes latitudes of grid. eg: c(-1,-5)
#' @param res The resolution of grid in degree
#' @param n.station Number of stations used per point for interpolation. The default is 8.
#' @param p Power parameter. The default is 1.
#'
#'
#' @return A data.frame with longitude, latitude and interpoled points
#'
#' @author Rodrigo Lins R. Jr., Fabricio  Daniel S. S. 
#'
#' @examples 
#' data(TempBrazil) # Temperature for some poins of Brazil
#'
#'
#' LonLat=TempBrazil[,1:2] #Data.frame with Longtude and Latitude
#' Temp=TempBrazil[,3] # Vector with observations in points
#'
#' LonInterval=c(-78,-34.10)  # Coordinates of extremes poins of longitude to grid
#' LatInterval=c(-36,5)  # Coordinates of extremes poins of latitude to grid
#'
#' Interpoled=IDW(xy=LonLat,z=Temp,xrange = LonInterval,yrange = LatInterval)
#' 
#' @export
IDW = function(xy,z,xrange,yrange,res=0.5,n.station=5,p=1){
  
  x.range <- seq(xrange[1],xrange[2],by=res) 
  y.range <- seq(yrange[1],yrange[2],by=res)
  
  coord=xy  
  INTERPOLED=data.frame()
  
  i=1
  j=1
  while (i<=length(x.range)) {
    
    while (j<=length(y.range)) {
      
      x1=x.range[i]          #Coordenada X do poonto a estimar
      y1=y.range[j]          #Coordenada y do ponto a estimar
      
      x2=coord[,1]           #Coordenadas x dos pontos observados e usados para estimar
      y2=coord[,2]           #Coordenadas y dos pontos observados e usados para estimar
      
      point=data.frame(x1,y1) # Data.frame com as coordenadas do ponto a estimar
      colnames(point)=c("x","y") 
      
      est_coordenadas=data.frame(x2,y2) #Data.frame com as coordenadas dos pontos observados
      colnames(est_coordenadas)=c("x","y")
      
      Di=distance(from = point,to=est_coordenadas) #Dist?ncia euclidianda entre os pontos
      
      aux=data.frame(coord,Di,z)
      colnames(aux)=c("Lat","Lon","Di","DATA")
      
      cron=aux[order(aux$Di),]
      cron=cron[1:n.station,]
      
      W=(1/(cron$Di**p))
      
      W[is.na(W)]=0
      
      PE=sum(W*cron$DATA)/sum(W)
      
      aux=data.frame(x1,y1,PE)
      INTERPOLED=rbind(aux,INTERPOLED)
      j=j+1
    }
    j=1
    i=i+1
  }
  
  colnames(INTERPOLED)=c("Longitude","Latitude","Variable")
  return(INTERPOLED)
}
