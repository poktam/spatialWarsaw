# dodać jakieś opisy
# testA
# sprawdzić czy testA oraz testB potrzebny poza pakietem (czy tylko funkcje wewnętrzne/pomocnicze) - ew. nazwy do zmiany
#' A TU BĘDZIE CAŁA POMOC DO FUNKCJI testA - chyba, że je ukryjemy i nie będą dostępne z zewnątrz ;)
#' @export
testA<-function(crds, clusters, sep, r_p=0.001, eps_r=10e-16, eps_np=10e-3, minPts0=16){
  params = as.list(environment())[-1] # zwraca listę parametrów funkcji bez pierwszego - czyli środowisko wewnątrz funkcji
  minPts=minPts0
  groups=testB(crds, minPts, sep, r_p, eps_r, eps_np)$cluster # run testB with initial parameters
  minPts_prev<-minPts
  if(length(unique(groups))>clusters+1) minPts=2*minPts else minPts=round(minPts/2,0) # adjust minPts
  groups_prev<-groups

  while(length(unique(groups))!=clusters+1 & minPts_prev!=minPts){ # repeat below until the desired number of clusters is obtained
    minPts=minPts+1

    groups=testB(crds, minPts, sep, r_p, eps_r, eps_np)$cluster # run testB

    if(length(unique(groups))>clusters+1 & length(unique(groups_prev))>clusters+1) { # adjust minPts
      minPts_prev=minPts
      minPts=2*minPts } else if(length(unique(groups))<clusters+1 & length(unique(groups_prev))<clusters+1){
        minPts_prev=minPts
        minPts=round(minPts/2,0)} else {
          rr=minPts
          minPts=round((minPts+minPts_prev)/2,0)
          minPts_prev=rr
        }
    groups_prev<-groups

  }
  message(paste0("minPts = ", minPts)) # Czy to zostawić?
  if(length(unique(groups))!=clusters+1) warning(paste0("Failed to get ", clusters, " clusters. Try changing minPts0 value."))

  # zastanowić się, czy klasa testing nie jest do usunięcia?
  structure(
    list(
      type = "A",
      cluster = groups,
      parameters = params
    ),
    class = "testing"
  )

}



# testB:
# sprawdzić czy testA oraz testB potrzebny poza pakietem (czy tylko funkcje wewnętrzne/pomocnicze) - ew. nazwy do zmiany
#' A TU BĘDZIE CAŁA POMOC DO FUNKCJI testB - chyba, że je ukryjemy i nie będą dostępne z zewnątrz ;)
#' @export
testB<-function(crds, minPts=5, sep=c(0.8, 0.6, 0.4, 0.2), r_p=0.001, eps_r=10e-16, eps_np=10e-3,
                bootstrap=FALSE, sample_size, times){
  params <- as.list(environment())[-1]
  output<-rep(length(sep)+1,nrow(crds))
  sep<-sort(sep, decreasing=TRUE)

  r<-r_p

  if(bootstrap){

    for(i in 1:length(sep)) {

      result0<-bootdbscan(crds, sample_size=sample_size, times=times, eps=r, minPts=minPts)$cluster # run bootdbscan with initial parameters
      noisePercent<-length(result0[result0==0])/nrow(crds)
      r_prev=r

      if(noisePercent>sep[i]) r=2*r else r=r/2 # adjust the radius
      nP_prev<-noisePercent

      while((abs(r-r_prev)>eps_r) & (abs(noisePercent-sep[i])>eps_np)){ # repeat below until the desired percentage of noise is obtained
        result0<-bootdbscan(crds, sample_size=sample_size, times=times, eps=r, minPts=minPts)$cluster # run bootdbscan
        noisePercent<-length(result0[result0==0])/nrow(crds)

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

      output<-ifelse(output==length(sep)+1 & result0>0, i, output) # final result
    }

  }else {

    for(i in 1:length(sep)) {
      result0<-dbscan(crds, r, minPts=minPts) # run dbscan
      noisePercent<-length(result0$cluster[result0$cluster==0])/nrow(crds)
      r_prev=r

      if(noisePercent>sep[i]) r=2*r else r=r/2 # adjust the radius
      nP_prev<-noisePercent

      while(abs(r-r_prev)>eps_r & abs(noisePercent-sep[i])>eps_np){ # repeat the following until the desired percentage of noise is obtained
        result0<-dbscan(crds, r, minPts=minPts) # run dbscan
        noisePercent<-length(result0$cluster[result0$cluster==0])/nrow(crds)

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

      output<-ifelse(output==length(sep)+1 & result0$cluster>0, i, output)} # adjust the radius

    if(length(sep)==1) output<-result0$cluster}

  # zastanowić się, czy klasa testing nie jest do usunięcia?
  structure(
    list(
      type = "B",
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
  cat("\nSummary of assignments:\n")
  print(summary(factor(x$cluster)))
  cat("\nAvailable fields: type, cluster, parameters\n")
}
