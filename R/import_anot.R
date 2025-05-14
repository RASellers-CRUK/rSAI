xml_to_list<-function(file) {
  result<-XML::xmlParse(file = file)
  rootnode<-XML::xmlRoot(result)
  xmlist<-XML::xmlToList(rootnode)
  return(xmlist)
}

#' Import an Inclusion Annotation
#'
#' This function loads a region annotation file produced in HALO and edits a Seurat object metadata,
#' displaying the overlap of every spot with the region/s and classifies each spot as within or 
#' outside of a region based upon a supplied cutoff. This function requires the spot annotation
#' produced by an "export_spot_annotation" family function. Spot annotations must match the supplied
#' Seurat object.
#'
#' @param seurat_data Seurat data to be annotated
#' @param spot_polygon Spot polygon annotations file (filepath) with multilayers produced by "export_spot_annotation" family function (rSAI)
#' @param region_polygon Region polygon annotations file (filepath) with multiple layers (if multiple regions) produced by HALO
#' @param overlap Percentage overlap required to include a spot within an annotation region (default=20)
#' @return The submitted Seurat object with annotated regions added to metadata
#' @export
import_include_region_Visium<-function(seurat_data,spot_polygon,region_polygon,overlap=20) {

  ## Reading in spot annotations
  xmlist<-xml_to_list(spot_polygon)
  coordlist<-list()
  for(annot in 1:length(xmlist)) {
    dat<-xmlist[[annot]]$Regions$Region
    coordlist[[xmlist[[annot]]$.attrs["Name"]]]<-as.data.frame(do.call(rbind,dat$Vertices))
  }
  if(any(is.na(match(colnames(seurat_data),names(coordlist))) | length(colnames(seurat_data))!=length(names(coordlist)))) {
    stop("Incompatible Seurat object and spot annotations; if filtering has been performed please rerun export_spot_annotation_fromSeurat")
  } else {
    ## Reading in inclusion annotation
    xmlist<-xml_to_list(region_polygon)
    regionlist<-list()
    for(annot in 1:length(xmlist)) {
      dat<-xmlist[[annot]]$Regions$Region
      regionlist[[gsub(" ","_",xmlist[[annot]]$.attrs["Name"])]]<-as.data.frame(do.call(rbind,dat$Vertices))
    }
    region_poly<-list()
    for(region in 1:length(regionlist)) {
      region_poly[[names(regionlist)[region]]]<-sp::SpatialPolygons(
        list(
          sp::Polygons(
            list(
              sp::Polygon(
                matrix(c(as.numeric(unlist(regionlist[[region]]["X"])),
                         as.numeric(unlist(regionlist[[region]]["Y"]))),ncol=2),hole=F)),ID = names(regionlist)[region])))
    }
    ## Calculate percentage overlap (each spot with each region)
    ol_list<-list()
    for(coord in 1:length(coordlist)) {
      coord_poly<-sp::SpatialPolygons(
        list(
          sp::Polygons(
            list(
              sp::Polygon(
                matrix(c(as.numeric(unlist(coordlist[[coord]]["X"])),
                         as.numeric(unlist(coordlist[[coord]]["Y"]))),ncol=2),hole=F)),ID = names(coordlist)[coord])))
      coord_ol<-list(ID=names(coordlist)[coord],
                     area=raster::area(coord_poly))
      for(region in 1:length(region_poly)) {
        inter<-raster::intersect(coord_poly,region_poly[[region]])
        if(is.null(inter)) {
          coord_ol[[names(region_poly)[region]]]<-0
        } else {
          coord_ol[[names(region_poly)[region]]]<-(raster::area(inter)/coord_ol[["area"]])*100
        }
      }
      ol_list[[names(coordlist)[coord]]]<-coord_ol
    }
    ol_df<-as.data.frame(do.call(rbind,ol_list))
    ol_df<-ol_df[match(colnames(seurat_data),rownames(ol_df)),]
    for(region in 1:length(regionlist)) {
      seurat_data[[paste0("pctOverlap_",names(regionlist)[region])]]<-unlist(ol_df[,names(regionlist)[region]])
      seurat_data[[paste0("pctOverlap",overlap,"_grp_",names(regionlist)[region])]]<-ifelse(ol_df[,names(regionlist)[region]]<=overlap,"out","in")
    }
    return(seurat_data)
  }
}

#' Import an Exclusion Annotation
#'
#' This function loads a region annotation file produced in HALO and subsets a Seurat object to
#' include only spots classfied as outside of the supplied annotation region/s based upon a supplied cutoff.
#' This function requires the spot annotation produced by an "export_spot_annotation" family function.
#' Spot annotations must match the supplied Seurat object.
#'
#' @param seurat_data Seurat data to be filtered
#' @param spot_polygon Spot polygon annotations file (filepath) with multilayers produced by rSAI
#' @param region_polygon Region polygon annotations file (filepath) with multiple layers (if multiple regions) produced by HALO
#' @param overlap Percentage overlap required to include a spot within an annotation region, marking the spot for filtration (default=70)
#' @return The submitted Seurat object filtered to remove spots within annotation regions
#' @export
import_exclude_region_Visium<-function(seurat_data,spot_polygon,region_polygon,overlap=70) {

  ## Reading in spot annotations
  xmlist<-xml_to_list(spot_polygon)
  coordlist<-list()
  for(annot in 1:length(xmlist)) {
    dat<-xmlist[[annot]]$Regions$Region
    coordlist[[xmlist[[annot]]$.attrs["Name"]]]<-as.data.frame(do.call(rbind,dat$Vertices))
  }
  if(any(is.na(match(colnames(seurat_data),names(coordlist))) | length(colnames(seurat_data))!=length(names(coordlist)))) {
    stop("Incompatible Seurat object and spot annotations; if filtering has been performed please rerun export_spot_annotation_fromSeurat")
  } else {
    ## Reading in exclusion annotation
    xmlist<-xml_to_list(region_polygon)
    excl<-list()
    for(annot in 1:length(xmlist)) {
      dat<-xmlist[[annot]]$Regions$Region
      excl[["exclude_region"]]<-as.data.frame(do.call(rbind,dat$Vertices))
    }
    excl_poly<-list()
    for(exc in 1:length(excl)) {
      excl_poly[[names(excl)[exc]]]<-sp::SpatialPolygons(
        list(
          sp::Polygons(
            list(
              sp::Polygon(
                matrix(c(as.numeric(unlist(excl[[exc]]["X"])),
                         as.numeric(unlist(excl[[exc]]["Y"]))),ncol=2),hole=F)),ID = names(excl)[exc])))
    }
    ol_list<-list()
    for(coord in 1:length(coordlist)) {
      coord_poly<-sp::SpatialPolygons(
        list(
          sp::Polygons(
            list(
              sp::Polygon(
                matrix(c(as.numeric(unlist(coordlist[[coord]]["X"])),
                         as.numeric(unlist(coordlist[[coord]]["Y"]))),ncol=2),hole=F)),ID = names(coordlist)[coord])))
      coord_ol<-list(ID=names(coordlist)[coord],
                     area=raster::area(coord_poly))
      for(exc in 1:length(excl_poly)) {
        inter<-raster::intersect(coord_poly,excl_poly[[exc]])
        if(is.null(inter)) {
          coord_ol[[names(excl_poly)[exc]]]<-0
        } else {
          coord_ol[[names(excl_poly)[exc]]]<-(raster::area(inter)/coord_ol[["area"]])*100
        }
      }
      ol_list[[names(coordlist)[coord]]]<-coord_ol
    }
    ol_df<-as.data.frame(do.call(rbind,ol_list))
    ol_df<-ol_df[match(colnames(seurat_data),rownames(ol_df)),]
    for(exc in 1:length(excl)) {
      excl_ol<-unlist(ol_df[,names(excl)[exc]])
      excl_ol_grp<-ifelse(ol_df[,names(excl)[exc]]>=overlap,"out","in")
    }
    seurat_data<-seurat_data[,colnames(seurat_data)[which(excl_ol_grp=="in")]]
    return(seurat_data)
  }
}
