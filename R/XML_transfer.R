xml_to_list<-function(file) {
  result<-XML::xmlParse(file = file)
  rootnode<-XML::xmlRoot(result)
  xmlist<-XML::xmlToList(rootnode)
  return(xmlist)
}

write_XML_layered<-function(coordlist,prefix,parallel=T) {
  hold<-c('<Annotations>')
  if(parallel) {
    anns<-parallel::mclapply(1:length(coordlist),function(x) {
      c(paste0('  <Annotation LineColor="65280" Name="',names(coordlist[x]),'" Visible="True">'),
        '    <Regions>',
        '      <Region Type="Polygon" HasEndcaps="0" NegativeROA="0">',
        '        <Vertices>',
        paste0('          <V X="',coordlist[[x]][,1],'" Y="',coordlist[[x]][,2],'" />'),
        '        </Vertices>',
        '        <Comments />',
        '      </Region>',
        '    </Regions>',
        '  </Annotation>')
    })
  } else {
    anns<-lapply(1:length(coordlist),function(x) {
      c(paste0('  <Annotation LineColor="65280" Name="',names(coordlist[x]),'" Visible="True">'),
        '    <Regions>',
        '      <Region Type="Polygon" HasEndcaps="0" NegativeROA="0">',
        '        <Vertices>',
        paste0('          <V X="',coordlist[[x]][,1],'" Y="',coordlist[[x]][,2],'" />'),
        '        </Vertices>',
        '        <Comments />',
        '      </Region>',
        '    </Regions>',
        '  </Annotation>')
    })
  }
  hold<-c(hold,
          unlist(anns),
          '</Annotations>')
  cat(hold,file = paste0(prefix,".annotations"),sep = "\n")
}

write_XML_1layer<-function(coordlist,prefix,parallel=T) {
  hold<-c('<Annotations>',
          '  <Annotation LineColor="65280" Name="Annotation Layer" Visible="True">',
          '    <Regions>')
  if(parallel) {
    anns<-parallel::mclapply(1:length(coordlist),function(x) {
      c('      <Region Type="Polygon" HasEndcaps="0" NegativeROA="0">',
        '        <Vertices>',
        paste0('          <V X="',coordlist[[x]][,1],'" Y="',coordlist[[x]][,2],'" />'),
        '        </Vertices>',
        '        <Comments />',
        '      </Region>')
    })
  } else {
    anns<-lapply(1:length(coordlist),function(x) {
      c('      <Region Type="Polygon" HasEndcaps="0" NegativeROA="0">',
        '        <Vertices>',
        paste0('          <V X="',coordlist[[x]][,1],'" Y="',coordlist[[x]][,2],'" />'),
        '        </Vertices>',
        '        <Comments />',
        '      </Region>')
    })
  }
  
  hold<-c(hold,
          unlist(anns),
          '    </Regions>',
          '  </Annotation>',
          '</Annotations>')
  cat(hold,file = paste0(prefix,"_onelayer.annotations"),sep = "\n")
}

geoJSON_op<-function(coordlist,outfile,parallel) {
  init<-geojsonR::TO_GeoJson$new()
  if(parallel) {
    flist<-parallel::mclapply(1:length(coordlist),function(x) {
      list(id=names(coordlist[x]),
           geometry=list(Polygon=list(split(as.matrix(coordlist[[x]]),seq(nrow(coordlist[[x]]))))),
           properties=list(isLocked="false",
                           measurements="[]",
                           classification=list(
                             name="Other",
                             colorRGB=-377282)
           )
      )
    })
  } else {
    flist<-lapply(1:length(coordlist),function(x) {
      list(id=names(coordlist[x]),
           geometry=list(Polygon=list(split(as.matrix(coordlist[[x]]),seq(nrow(coordlist[[x]]))))),
           properties=list(isLocked="false",
                           measurements="[]",
                           classification=list(
                             name="Other",
                             colorRGB=-377282)
           )
      )
    })
  }
  names(flist)<-rep("Feature",length(flist))
  fdat<-list(features=flist)
  fobj<-init$FeatureCollection(fdat, stringify = TRUE)
  
  cat(fobj$json_dump,file = outfile)
}

expoi_oline<-function(bcs,tissue,radius) {
  
  XY<-data.frame(Var1=bcs$imagecol[tissue],Var2=bcs$imagerow[tissue])
  subcoords<-expand.grid((XY$Var1-ceiling(radius)):(XY$Var1+ceiling(radius)),
                         (XY$Var2-ceiling(radius)):(XY$Var2+ceiling(radius)))
  
  cdist<-raster::pointDistance(XY,subcoords,lonlat = F)
  poi<-subcoords[which(cdist<=radius),]
  poi2<-subcoords[which(round(cdist)==round(radius)),]
  
  hold<-c()
  hold2<-c(rownames(poi2)[1])
  for(x in 1:nrow(poi2)) {
    dists<-raster::pointDistance(poi2[hold2[x],],poi2,lonlat=F)
    hold<-append(hold,which.min(dists),after=length(hold))
    hold2<-append(hold2,names(which.min(dists[!names(dists) %in% names(hold)])),after=length(hold2))
  }
  poi2<-poi2[hold2,]
  oline<-rbind(poi2,poi2[1,])
  colnames(oline)<-c("x","y")
  
  return(list(sample=bcs$barcode[tissue],
              oline=oline))
}

expoi_corners<-function(bcs,tissue,radius) {
  
  XY<-data.frame(x=bcs$imagecol[tissue],y=bcs$imagerow[tissue])
  TL<-c(round(XY[,"x"]-radius),round(XY[,"y"]+radius))
  BR<-c(round(XY[,"x"]+radius),round(XY[,"y"]-radius))
  corners<-rbind(TL,BR)
  colnames(corners)<-c("x","y")
  
  return(list(sample=bcs$barcode[tissue],
              corners=corners))
}

geoJSON_opL<-function(coordlist,outfile,parallel) {
  init<-geojsonR::TO_GeoJson$new()
  if(parallel) {
    flist<-parallel::mclapply(1:length(coordlist),function(x) {
      list(id=names(coordlist[x]),
           geometry=list(Polygon=list(split(as.matrix(coordlist[[x]]),seq(nrow(coordlist[[x]]))))),
           properties=list(isLocked="false",
                           measurements="[]",
                           classification=list(
                             name="Other",
                             colorRGB=-377282)
           )
      )
    })
  } else {
    flist<-lapply(1:length(coordlist),function(x) {
      list(id=names(coordlist[x]),
           geometry=list(Polygon=list(split(as.matrix(coordlist[[x]]),seq(nrow(coordlist[[x]]))))),
           properties=list(isLocked="false",
                           measurements="[]",
                           classification=list(
                             name="Other",
                             colorRGB=-377282)
           )
      )
    })
  }
  names(flist)<-rep("Feature",length(flist))
  fdat<-list(features=flist)
  
  fobjs<-lapply(flist,function(obj) {init$Feature(obj, stringify = TRUE)})
  
  cat("[",file = outfile)
  for(i in 1:(length(fobjs)-1)) {
    cat(paste0(fobjs[[i]]$json_dump,","),file = outfile,append = T)
  }
  cat(fobjs[[length(fobjs)]]$json_dump,file = outfile,append = T)
  cat("]",file = outfile,append = T)
}


#' Convert 1layer Annotations to Multilayer
#'
#' This function loads a spot annotation file in 1layer format and uses the spot names
#' from a submitted Seurat data object to convert the spot annotations to multilayer format.
#' This function is intended for use once annotations have been registered using HALO from
#' the 1layer spot annotations but can also be used to convert any 1layer format spot
#' annotations file in to multilayer format with matched spots.
#'
#' @param seurat_data Seurat data to be annotated
#' @param outfile Name of output file
#' @param spot_anot Spot polygon annotations file (filepath) with multilayers produced by "export_spot_annotation" family function (rSAI)
#' @param reg_spot_anot 1layer format spot annotations file (filepath) to convert to multilayer format
#' @param parallel Use parallelisation to speed up function (recommended for multi-core processors)
#' @export
convert_registered_1layer<-function(seurat_data,spot_anot,reg_spot_anot,parallel=T) {
  prefix<-sub(".annotations$","_converted",reg_spot_anot)
  ## Reading in spot annotations
  xmlist<-xml_to_list(spot_anot)
  coordlist<-list()
  if(length(xmlist)==ncol(seurat_data)) {
    for(annot in 1:length(xmlist)) {
      dat<-xmlist[[annot]]
      coordlist[[colnames(seurat_data)[annot]]]<-as.data.frame(do.call(rbind,dat$Regions$Region$Vertices))
    }
  } else {
    stop("Supplied Seurat data not matching spot annotations")
  }
  ## Reading in registered spot annotations (1layer)
  xmlist<-xml_to_list(reg_spot_anot)
  vlist<-list()
  if(length(xmlist$Annotation[[1]])==ncol(seurat_data)) {
    for(annot in 1:length(xmlist$Annotation[[1]])) {
      dat<-xmlist$Annotation[[1]][[annot]]
      vlist[[names(coordlist)[annot]]]<-as.data.frame(do.call(rbind,dat$Vertices))
    }
  } else {
    stop("Registered spot annotations not matching spot annotations")
  }
  write_XML_layered(vlist,prefix,parallel)
}

#' Split Annotations File in to Layers
#'
#' This function takes a table of barcodes with their allocated regions along with a spot
#' annotation file and produces a filtered spot annotation file containing only spots
#' within the defined regions. Spots within the same region are stored within 1 layer.
#' This function is intended for visualisation of regions within HALO.
#'
#' @param spot_anot Spot polygon annotations file (filepath) with multilayers produced by "export_spot_annotation" family function (rSAI)
#' @param region_df Region data.frame containing two columns: col1 - barcodes/spot names, col2 - allocated region (no region = NA)
#' @param outfile Name of output file
#' @export
split_annotations<-function(spot_anot,region_df,outfile) {
  
  ## Reading in spot annotations
  xmlist<-xml_to_list(spot_anot)
  vlist<-list()
  for(annot in 1:length(xmlist)) {
    dat<-xmlist[[annot]]
    vlist[[dat$.attrs["Name"]]]<-as.data.frame(do.call(rbind,dat$Regions$Region$Vertices))
  }
  if(length(which(is.na(match(unique(region_df$barcode),names(vlist)))))>0) {
    stop("Unmatching barcodes supplied in region data.frame")
  } else {
    colnames(region_df)<-c("barcode","region")
    region_df<-region_df[order(region_df$region),]
    region_list<-split(region_df,region_df$region)
    ## Output selected grad spots as region layers
    cat("<Annotations>\n",file = outfile)
    for(reg in 1:length(region_list)) {
      cat('  <Annotation LineColor="65280" Name="',names(region_list)[reg],'" Visible="True">\n',sep = "",file = outfile,append = T)
      cat('    <Regions>\n',sep = "",file = outfile,append = T)
      for(x in 1:nrow(region_list[[reg]])) {
        cat('      <Region Type="Polygon" HasEndcaps="0" NegativeROA="0">\n',sep = "",file = outfile,append = T)
        cat('        <Vertices>\n',file = outfile,append = T)
        bcode<-region_list[[reg]]$barcode[x]
        verts<-vlist[[bcode]]
        for(i in 1:nrow(verts)) {
          cat('          <V X="',verts[i,"X"],'" Y="',verts[i,"Y"],'" />\n',sep = "",file = outfile,append = T)
        }
        cat('        </Vertices>\n',sep = "",file = outfile,append = T)
        cat('        <Comments />\n',sep = "",file = outfile,append = T)
        cat('      </Region>\n',sep = "",file = outfile,append = T)
      }
      cat('    </Regions>\n',sep = "",file = outfile,append = T)
      cat('  </Annotation>\n',sep = "",file = outfile,append = T)
    }
    cat("</Annotations>",sep = "",file = outfile,append = T)
  }
}

#' Produce Spot Annotation Polygons from Space Ranger Outputs
#'
#' This function provides a method for producing HALO-compatible annotation files
#' containing polygon vertices (and/or ellipse vertices) for 10x Visium spot data.
#' The inputs for this function are the Space Ranger outputs along with the
#' full-resolution TIFF image used as input to Space Ranger. The TIFF image is not
#' loaded in to memory and is only used for extracting the tagged image resolution
#' this is necessary as the defined spot radii defined by Space Ranger are inaccurate
#' leading to spot diameters of 65um (where 55um is the manufacturer definition; current
#' for Space Ranger version 1.3). Annotations using multilayer format are produced where
#' each spot resides within a named layer of the annotation file; 1layer format spot
#' annotations are intended for use when registering coordinates between serial-sections
#' (the 1layer coordinates are registered using HALO and then reconverted to multilayer
#' format via "convert_registered_1layer"). This function may take a while to run but
#' should only need to be ran once per sample.
#'
#' @param folderpath Folderpath of Space Ranger outputs (path)
#' @param imagepath Filepath of full-resolution TIFF image used for Space Ranger input
#' @param prefix A prefix used for the output annotation files; can include a directory structure
#' @param ellipse Output annotations using ellipse shape format rather than polygon format; faster but not compatible with region functions
#' @param onelayer Output both multilayer and 1layer format annotation files; needed for serial-section registration
#' @param parallel Use parallelisation to speed up function (recommended for multi-core processors)
#' @export
export_spot_annotation_fromBase<-function(folderpath,imagepath,prefix="sample",ellipse=F,onelayer=T,parallel=T) {

  alldirs<-list.dirs(folderpath,recursive = T)
  sdirs<-dirname(alldirs[grep("spatial$",alldirs,perl=T)])
  
  tifftags<-suppressWarnings(tiff::readTIFF(imagepath,payload=F))
  tl<-XML::xmlToList(tifftags$description)
  resol<-mean(as.numeric(c(tl$Image$Pixels$.attrs["PhysicalSizeX"],tl$Image$Pixels$.attrs["PhysicalSizeY"])))
  diameter<-55/resol
  radius<-diameter/2
  
  scalefactor_paths<-grep("scalefactors_json",list.files(paste(sdirs,"spatial",sep="/"),full.names = T),value=T)
  tissue_paths<-grep("tissue_positions",list.files(paste(sdirs,"spatial",sep="/"),full.names = T),value=T)
  
  scales<-rjson::fromJSON(file = scalefactor_paths)
  bcs<-read.csv(tissue_paths,col.names=c("barcode","tissue","row","col","imagerow","imagecol"), header = FALSE)
  bcs<-bcs[order(bcs$barcode),]
  bcs$tissue<-as.factor(bcs$tissue)
  bcs<-bcs[which(bcs$tissue==1),]
  bcs<-data.frame(barcode=bcs$barcode,
                  tissue=bcs$tissue,
                  sapply(bcs[,c("row","col","imagerow","imagecol")],as.numeric))
  
  if(!ellipse) {
    if(parallel) {
      res<-parallel::mclapply(which(bcs$tissue==1),function(a){expoi_oline(bcs,a,radius)})
    } else {
      res<-lapply(which(bcs$tissue==1),function(a){expoi_oline(bcs,a,radius)})
    }
    coordlist<-list()
    for(x in 1:length(res)) {
      coordlist[[res[[x]]$sample]]<-res[[x]]$oline
    }
    rm(res)
    saveRDS(coordlist,paste0(prefix,"_coordlist.rds"))
    write_XML_layered(coordlist,prefix,parallel)
    if(onelayer) {
      write_XML_1layer(coordlist,prefix,parallel)
    }
  } else {
    if(parallel) {
      res<-parallel::mclapply(which(bcs$tissue==1),function(a){expoi_corners(bcs,a,radius)})
    } else {
      res<-lapply(which(bcs$tissue==1),function(a){expoi_corners(bcs,a,radius)})
    }
    coordlist<-list()
    for(x in 1:length(res)) {
      coordlist[[res[[x]]$sample]]<-res[[x]]$corners
    }
    rm(res)
    saveRDS(coordlist,paste0(prefix,"_coordlist_ellipse.rds"))
    write_XML_layered(coordlist,prefix,parallel)
    if(onelayer) {
      write_XML_1layer(coordlist,prefix,parallel)
    }
  }
}

#' Reproduce Spot Annotations using Filtered Seurat Object
#'
#' This function is intended for use once a base spot annotations file has been produced
#' by "export_spot_annotation_fromBase" and any filtration of spots has occurred. This
#' function takes as input the base spot annotations along with the filtered Seurat object
#' and outputs a filtered multilayer spot annotation file (allowing also the output of
#' 1layer format annotations).
#'
#' @param seurat_data Filtered Seurat object to be used for the filtered output
#' @param base_spot_anot Base spot annotation file (multilayer) produced by "export_spot_annotation_fromBase"
#' @param image_name Name of the image coordinates within the Seurat object to be used (default is "slice1"; `seurat_data@images[["slice1"]]`)
#' @param prefix A prefix used for the output annotation files; can include a directory structure
#' @param onelayer Output both multilayer and 1layer format annotation files; needed for serial-section registration
#' @param parallel Use parallelisation to speed up function (recommended for multi-core processors)
#' @export
export_spot_annotation_fromSeurat<-function(seurat_data,base_spot_anot,image_name="slice1",prefix="sample_filtered",onelayer=T,parallel=T) {
  
  xmlist<-xml_to_list(base_spot_anot)
  coordlist<-list()
  for(annot in 1:length(xmlist)) {
    dat<-xmlist[[annot]]$Regions$Region
    coordlist[[xmlist[[annot]]$.attrs["Name"]]]<-as.data.frame(do.call(rbind,dat$Vertices))
  }
  
  bcs<-seurat_data@images[[image_name]]@coordinates
  bcs$tissue<-as.factor(bcs$tissue)
  bcs<-bcs[which(bcs$tissue==1),]
  
  if(length(which(is.na(match(rownames(bcs),names(coordlist)))))>0) {
    stop("Unmatched barcodes in seurat_data")
  } else {
    
    coordlist<-coordlist[rownames(bcs)]
    write_XML_layered(coordlist,prefix,parallel)
    if(onelayer) {
      write_XML_1layer(coordlist,prefix,parallel)
    }
  }
}

#' Convert Annotations to QuPath geoJSON file
#'
#' This function converts a HALO-compatible .annotations file to a QuPath-compatible geoJSON file.
#'
#' @param annotations HALO-compatible .annotations file either created in HALO or by rSAI
#' @param outfile Output file name
#' @param toJSONlist File format (list of JSON features [default/T] OR FeatureCollection)
#' @param parallel Use parallelisation to speed up function (recommended for multi-core processors)
#' @export
convert_annot_to_geoJSON<-function(annotations,outfile,toJSONlist=T,parallel=T) {
  
  ## Read in annotations
  xmlist<-xml_to_list(annotations)
  coordlist<-list()
  for(annot in 1:length(xmlist)) {
    dat<-xmlist[[annot]]
    coordlist[[xmlist[[annot]]$.attrs["Name"]]]<-sapply(as.data.frame(do.call(rbind,dat$Regions$Region$Vertices)),as.numeric)
  }
  ## Output geoJSON
  if(toJSONlist) {
    geoJSON_opL(coordlist,outfile,parallel)
  } else {
    geoJSON_op(coordlist,outfile,parallel) 
  }
}

#' Convert QuPath geoJSON file to Annotations
#'
#' This function converts a QuPath-compatible geoJSON file to a HALO-compatible .annotations file.
#'
#' @param annotations QuPath-compatible geoJSON file
#' @param outfile Output file name
#' @param fromJSONlist File format (list of JSON features [default/T] OR FeatureCollection)
#' @param onelayer Output both multilayer and 1layer format annotation files; needed for serial-section registration
#' @param parallel Use parallelisation to speed up function (recommended for multi-core processors)
#' @export
convert_geoJSON_to_annot<-function(geoJSON,prefix="sample_json",fromJSONlist=T,onelayer=T,parallel=T) {
  
  if(fromJSONlist) {
    ## Read in geoJSON_fromJSONlist
    coord_json<-jsonlite::fromJSON(geoJSON)
    
    ## Wrangle to coordlist format [named list, each barcode has coords as dataframe with c("x","y") colnames]
    coordlist<-sapply(coord_json,function(z) {
      poly<-matrix(unlist(z$geometry$coordinates),ncol=2,byrow=T)
      data.frame(x=poly[,1],
                 y=poly[,2])
    },simplify = F)
    names(coordlist)<-sapply(coord_json,function(z) {
      z$id
    })
  } else {
    ## Read in geoJSON_fromFeatureCollectionJSON
    coord_json<-geojsonR::FROM_GeoJson(geoJSON)
    
    ## Wrangle to coordlist format [named list, each barcode has coords as dataframe with c("x","y") colnames]
    coordlist<-sapply(coord_json[[1]],function(z) {
      data.frame(x=z$geometry$coordinates[,1],
                 y=z$geometry$coordinates[,2])
    },simplify = F)
    names(coordlist)<-sapply(coord_json[[1]],function(z) {
      z$id
    })
  }
  write_XML_layered(coordlist,prefix,parallel)
  if(onelayer) {
    write_XML_1layer(coordlist,prefix,parallel)
  }
}

