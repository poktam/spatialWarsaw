# sprawdzić czy ClustConti() oraz ClustDisjoint() potrzebny poza pakietem
# (czy tylko funkcje wewnętrzne/pomocnicze) - ew. nazwy do zmiany
#
####################
### ClustConti() ###
####################
#

#' Linijka nr 1 - function title
#'
#' Linijka nr 2 - description
#'
#' Linijka nr 3 - A TU BĘDZIE CAŁA POMOC DO FUNKCJI ClustConti() - chyba, że je ukryjemy i nie będą dostępne z zewnątrz ;)
#'
#'
#' @name ClustConti
#' @param data_sf (obiekt sf lub data.frame zawierający współrzędne - w data frame 1 kolumna musi być X coords, druga kolumna Y coords)
#' @param clusters do opisu
#' @param sep do opisu
#' @param r_p do opisu
#' @param eps_r do opisu
#' @param eps_np do opisu
#' @param minPts0 do opisu
#'
#' @examples #To be done!!!
#'
#' @return `ClustConti()` returns ... to be done.
#'
#' @export
ClustConti<-function(data_sf, clusters, sep, r_p=0.001, eps_r=10e-16, eps_np=10e-3, minPts0=16){
  params <- as.list(environment())[-1]

  if((inherits(data_sf,"sf") && st_geometry_type(data_sf,FALSE)=="POINT")) {
    crds<-as.data.frame(st_coordinates(data_sf))
    colnames(crds)<-c("X_coord","Y_coord")
    cat("Data_sf was detected as an object of class sf.\n", sep = "")
  }  else if(inherits(data_sf,"data.frame",TRUE)==1){
    crds<-data_sf[,c(1,2)]
    colnames(crds)<-c("X_coord","Y_coord")
    cat("Data_sf was detected as an object of class data.frame.\n", sep = "")
  }  else {
    stop("The class of data_sf must only be 'sf' of geometry type 'POINTS' or 'data.frame'.")
  }

  # zbadać clusters, czy nie ujemne - może inny warunek niż <=0
  clusters<-round(clusters,0)
  if (clusters<=0) {
    clusters<-4
    cat("Incorrect clusters parameter. Set to the proposed value: 4.","\n",sep="")
  }

  # zbadać minPts0, czy nie ujemne - może inny warunek niż <=0
  minPts0<-round(minPts0,0)
  if (minPts0<=0) {
    minPts0<-16
    cat("Incorrect minPts0 parameter. Set to the proposed value: 16.","\n",sep="")
  }

  minPts<-minPts0
  sink(nullfile()) #ew. spróbować z capture.output() i invisible()
  output<-ClustDisjoint(data_sf=crds, sep=sep, r_p=r_p, eps_r=eps_r, eps_np=eps_np, minPts=minPts)$cluster # run ClustDisjoint with initial parameters
  sink()
  minPts_prev<-minPts

  if(length(unique(output))>clusters+1) {
    minPts=2*minPts
  } else {
    minPts=round(minPts/2,0)} # adjust minPts
  output_prev<-output

  while(length(unique(output))!=clusters+1 & minPts_prev!=minPts){ # repeat below until the desired number of clusters is obtained
    minPts=minPts+1

    sink(nullfile())
    output<-ClustDisjoint(data_sf=crds, sep=sep, r_p=r_p, eps_r=eps_r, eps_np=eps_np, minPts=minPts)$cluster # run ClustDisjoint
    sink()

    if(length(unique(output))>clusters+1 & length(unique(output_prev))>clusters+1) {
      minPts_prev=minPts
      minPts=2*minPts } else if(length(unique(output))<clusters+1 & length(unique(output_prev))<clusters+1){
        minPts_prev=minPts
        minPts=round(minPts/2,0)} else {
          rr=minPts
          minPts=round((minPts+minPts_prev)/2,0)
          minPts_prev=rr  # adjust minPts
        }
    output_prev<-output

  }
  cat("Final value: minPts = ", minPts,"\n",sep="")
  if(length(unique(output))!=clusters+1) warning(paste0("Failed to get ", clusters, " clusters. Try changing minPts0 value."))

  # zastanowić się, czy klasa testing nie jest do usunięcia?
  structure(
    list(
      type = "ClustConti",
      cluster = output,
      parameters = params
    ),
    class = "testing"
  )

}


#######################
### ClustDisjoint() ###
#######################
#
# sprawdzić czy ClustConti() oraz ClustDisjoint() potrzebny poza pakietem (czy tylko funkcje wewnętrzne/pomocnicze) - ew. nazwy do zmiany

#' Linijka nr 1 - function title
#'
#' Linijka nr 2 - description
#'
#' Linijka nr 3 - A TU BĘDZIE CAŁA POMOC DO FUNKCJI ClustDisjoint() - chyba, że je ukryjemy i nie będą dostępne z zewnątrz ;)
#'
#'
#' @name ClustDisjoint
#' @param data_sf (obiekt sf lub data.frame zawierający współrzędne - w data frame 1 kolumna musi być X coords, druga kolumna Y coords)
#' @param sep do opisu
#' @param r_p do opisu
#' @param eps_r do opisu
#' @param eps_np do opisu
#' @param minPts do opisu
#' @param bootstrap do opisu
#' @param sample_size do opisu (used only if bootstrap=TRUE)
#' @param times do opisu  (used only if bootstrap=TRUE)
#'
#' @examples #To be done!!!
#'
#' @return `ClustDisjoint()` returns ... to be done.
#'
#' @export
ClustDisjoint<-function(data_sf, sep=c(0.8, 0.6, 0.4, 0.2), r_p=0.001, eps_r=10e-16, eps_np=10e-3, minPts=5,
                bootstrap=FALSE, sample_size=NULL, times=NULL){
  params <- (as.list(environment()))[-1] # tu można poprawić, żeby w zależności od bootstrap dawało odpowiednie

  if((inherits(data_sf,"sf") && st_geometry_type(data_sf,FALSE)=="POINT")) {
    crds<-as.data.frame(st_coordinates(data_sf))
    colnames(crds)<-c("X_coord","Y_coord")
    cat("Data_sf was detected as an object of class sf.\n", sep = "")
  }  else if(inherits(data_sf,"data.frame",TRUE)==1){
    crds<-data_sf[,c(1,2)]
    colnames(crds)<-c("X_coord","Y_coord")
    cat("Data_sf was detected as an object of class data.frame.\n", sep = "")
  }  else {
    stop("The class of data_sf must only be 'sf' of geometry type 'POINTS' or 'data.frame'.")
  }

  # zbadać minPts, czy nie ujemne - może inny warunek niż <=0
  minPts<-round(minPts,0)
  if (minPts<=0) {
    minPts<-5
    cat("Incorrect minPts parameter. Set to the proposed value: 5.","\n",sep="")
  }

  if(bootstrap){

    if (is.null(sample_size) || is.null(times)) {
      stop("For bootstrap=TRUE, both sample_size and times parameters must be set.")
    }

    # zbadać sample_size i ustawić ew. na wielkość zbioru danych
    sample_size <- as.integer(sample_size)
    if (sample_size > nrow(crds) || sample_size<1) {
      sample_size<-nrow(crds)
      cat("Wrong sample size. Sample_size set to:",sample_size,"\n",sep="")
    }

    # zbadać times, czy nie ujemne - może inny warunek niż <=0
    times<-round(times,0)
    if (times<=0) {
      times<-16
      cat("Incorrect times value. Times set to the proposed value: 16.","\n",sep="")
    }
  }

  sep<-sort(sep, decreasing=TRUE)
  output<-rep(length(sep)+1,nrow(crds))

  r<-r_p

  # można będzie spróbować zmniejszyć kod w dużym if (poniżej), bo się powtarza
  if(bootstrap){
    for(i in 1:length(sep)) {
      sink(nullfile())
      result_temp<-bootdbscan(data_sf=crds, sample_size=sample_size, times=times, eps=r, minPts=minPts, plot=FALSE)$cluster # run bootdbscan with initial parameters
      sink()
      noisePercent<-length(result_temp[result_temp==0])/nrow(crds)
      r_prev=r

      if(noisePercent>sep[i]) r=2*r else r=r/2 # adjust the radius
      nP_prev<-noisePercent

      while((abs(r-r_prev)>eps_r) & (abs(noisePercent-sep[i])>eps_np)){ # repeat below until the desired percentage of noise is obtained
        sink(nullfile())
        result_temp<-bootdbscan(data_sf=crds, sample_size=sample_size, times=times, eps=r, minPts=minPts, plot=FALSE)$cluster # run bootdbscan
        sink()
        noisePercent<-length(result_temp[result_temp==0])/nrow(crds)

        if(noisePercent>sep[i] & nP_prev>sep[i]) { # adjust theradius
          r_prev=r
          r=2*r } else if(noisePercent<sep[i] & nP_prev<sep[i]){
            r_prev=r
            r=r/2} else {
              rr=r
              r=(r+r_prev)/2
              r_prev=rr}

        nP_prev<-noisePercent

      }

      output<-ifelse(output==length(sep)+1 & result_temp>0, i, output) # final result
    }

  }else {

    for(i in 1:length(sep)) {
      result_temp<-dbscan(crds, r, minPts=minPts) # run dbscan
      noisePercent<-length(result_temp$cluster[result_temp$cluster==0])/nrow(crds)
      r_prev=r

      if(noisePercent>sep[i]) r=2*r else r=r/2 # adjust the radius
      nP_prev<-noisePercent

      while(abs(r-r_prev)>eps_r & abs(noisePercent-sep[i])>eps_np){ # repeat the following until the desired percentage of noise is obtained
        result_temp<-dbscan(crds, r, minPts=minPts) # run dbscan
        noisePercent<-length(result_temp$cluster[result_temp$cluster==0])/nrow(crds)

        if(noisePercent>sep[i] & nP_prev>sep[i]) { # adjust the radius
          r_prev=r
          r=2*r } else if(noisePercent<sep[i] & nP_prev<sep[i]){
            r_prev=r
            r=r/2} else {
              rr=r
              r=(r+r_prev)/2
              r_prev=rr}

        nP_prev<-noisePercent
      }

      output<-ifelse(output==length(sep)+1 & result_temp$cluster>0, i, output)} # adjust the radius

    if(length(sep)==1) output<-result_temp$cluster}

  # zastanowić się, czy klasa testing nie jest do usunięcia?
  structure(
    list(
      type = "ClustDisjoint",
      cluster = output,
      parameters = params
    ),
    class = "testing"
  )

}

# Define the output format for print(testing)
# MK: Może do usunięcia jak testA i testB tylko wewnętrzne?
#' @export
print.testing <- function(x) {
  cat("Type", x$type,"group assignments for", length(x$cluster), "objects.\n")
  cat("Summary of assignments:\n")
  print(summary(factor(x$cluster)))
  cat("\nAvailable fields: type, cluster, parameters\n")
}
