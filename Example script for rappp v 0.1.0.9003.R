#### EXAMPLE SCRIPT USING rappp_0.1.0.9003
{ # Run from here for full script
  { # All loading and matrix naming, incl packages
    ################### LOAD Packages ############################################
    {
      # ------------------------ Packages & functions ---------------------------------
      library(rappp) # all ap_-functions
      library(SBA) # read_FlexMAP3D()
      library(WriteXLS) # WriteXLS()
      #-----------------------------------------------------------------------
    }

    ################### LOAD DATA & INFO #########################################
    { # Load and format all input
      # ------------------------ Basic Project input ---------------------------------
      {
        analysis <- "DRNA_" # I use as prefix for all filenames.
      }
      #-----------------------------------------------------------------------

      # -------- Luminex files ------------------------------------------
      {
        assay_1_org <- read_FlexMAP3D("Luminex files/DRNA_plate1_shipment2_CH_20190718_140301.csv")
        mfi_1 <- assay_1_org$median
        count_1 <- assay_1_org$count

        assay_2_org <- read_FlexMAP3D("Luminex files/DRNA_plate2_shipment2_CH_20190718_140820.csv")
        mfi_2 <- assay_2_org$median
        count_2 <- assay_2_org$count

        mfi <- rbind(mfi_1, mfi_2)
        count <- rbind(count_1, count_2)

        mfi_ct <- read.delim("Luminex files/AY144_data_intensity.tsv")
        mfi_ct <- mfi_ct[grep("^His", mfi_ct$sample.antigen), -1]
      }
      #-----------------------------------------------------------------------

      # -------- Information files ------------------------------------------
      {
        beads_org <- read.delim(paste0("Information files/AY144_antigen.tsv"), stringsAsFactors=F)

        sampleinfo_org <- read.delim("Information files/DRNA_01 samplelayout.txt", na.strings="", stringsAsFactors=F)
        sampleinfo <- sampleinfo_org
      }
      #-----------------------------------------------------------------------

      # -------- Adjusting matrices ------------------------------------------
      # Fix with gene-names in beads
      # PrESTs mapping to more than one gene are shortened to only the first gene followed by an asterix (*)
      # The new gene-names are combined with HPRR to use as an identifier for the rest of the script.
      {
        beads <- data.frame(beads_org,
                            GeneShort=unlist(lapply(strsplit(as.character(beads_org$Gene.name), ","),
                                                    function(y) ifelse(length(y) > 1, paste0(y[1],"*"),y[1]))), stringsAsFactors=F)
        beads <- data.frame(Gene_HPRR=ifelse(grepl("HPRA", beads$Antigen.name),
                                             paste0(beads$GeneShort,"_HPRR",beads$PrEST.ID),
                                             paste(beads$Antigen.name)),
                            beads, stringsAsFactors=F)
        rownames(beads) <- beads$Gene_HPRR
      }


      # Rename row- and colnames
      {
        rownames(sampleinfo) <- ifelse(grepl("EMPTY|MIX",sampleinfo$sample_name),
                                       as.character(sampleinfo$tube_label), # Replicates and buffer wells with conescutive numbering to allow as rownames.
                                       as.character(sampleinfo$sample_name))

        colnames(mfi) <- beads$Gene_HPRR ; rownames(mfi) <- rownames(sampleinfo)
        colnames(count) <- beads$Gene_HPRR ; rownames(count) <- rownames(sampleinfo)


        colnames(mfi_ct) <- beads$Gene_HPRR
      }
    } # END Loading and matrix naming
  } # END Loading and matrix naming, incl packages

  { # Create SBA-list and do QC
    #--------------- Create data lists --------------------------------------------
    #
    ## Combine data sets
    # The below structure is used by all used ap_-functions,
    # and is built upon if function output is stored in a new/the same variable.
    SBA <- list(MFI=mfi,
                COUNT=count,
                SAMPLES=sampleinfo,
                BEADS=beads,
                CT=mfi_ct,
                FILTERINFO=NULL)
    #
    #-------------------------------------------------------------------------

    # -------- QC plots and filtering ------------------------------------------
    #
    shouldplot=T
    ## Coupling efficiency test
    SBA <- ap_ct(SBA, grep("Bare", SBA$BEADS$Gene_HPRR), empty_co_multiple=3, shouldplot=shouldplot,
                 filename=paste0("Results/Second shipment/QC/",analysis,"Coupling-efficiency_",format(Sys.time(),"%Y-%m-%d_%H%M%S"),".pdf"),
                 width=30, height=6, useDingbats=F)

    ## Anti-human IgG
    SBA <- ap_igx(SBA, grep("IgG", SBA$BEADS$Gene_HPRR), IgType="G", IgX_cutoff=5000, cosfac=c(3, -3),
                  shouldplot=shouldplot,
                  filename=paste0("Results/Second shipment/QC/",analysis,"ahIgG_",format(Sys.time(),"%Y-%m-%d_%H%M%S"),".pdf"),
                  width=10, height=6, useDingbats=F)

    ## Bead count
    SBA <- ap_count(SBA, labels="Gene_HPRR", protein="GeneShort", agID="PrEST.ID",
                    samp_co=32, bead_flag=32, bead_filter=16, N_filter=0, shouldplot=shouldplot,
                    filename=paste0("Results/Second shipment/QC/",analysis,"BeadCount_",format(Sys.time(),"%Y-%m-%d_%H%M%S"),".pdf"),
                    width=12, height=10, useDingbats=F)

    ## Manual annotations
    SBA$SAMPLES$Filtered[grep("rep|pool|commercial|mix|blank|empty|buffer", SBA$SAMPLES$Sample, ignore.case=T)] <- "Control"
    SBA$BEADS$Filtered[grep("ebna|ahig|human", SBA$BEADS$Gene_HPRR, ignore.case=T)] <- "PosControl"
    SBA$BEADS$Filtered[grep("bare|empty|his6abp|hisabp", SBA$BEADS$Gene_HPRR, ignore.case=T)] <- "NegControl"

    # Summary of filtered and flaged
    print(table(SBA$BEADS$Filtered))
    print(table(SBA$BEADS$Flagged))
    print(table(SBA$SAMPLES$Filtered))
    #
    #-------------------------------------------------------------------------
  } # END QC filtering

  # -------- Data overview ------------------------------------------
  #
  {
    ## Signal overview
    ap_overview(SBA, filename=paste0("Results/Second shipment/QC/",analysis,"Signal-overview_",format(Sys.time(),"%Y-%m-%d_%H%M%S"),".pdf"),
                width=20, height=15, useDingbats=F, cex.axis=0.3)

    ## CV and replicates
    ap_rep(SBA, iter=500, filename=paste0("Results/Second shipment/QC/",analysis,"reproducibility_",format(Sys.time(),"%Y-%m-%d_%H%M%S"),".pdf"),
           width=12, height=12, useDingbats=F)

    ## tSNE
    tsne_perp(list(SBA$MFI), perp=c(2,5,10,50), iterations=1000,
              groups=as.factor(SBA$SAMPLES$AssayNum),
              names=rownames(SBA$MFI),
              legend=T,
              filename=paste0("Results/Second shipment/QC/",analysis,"tSNE_",format(Sys.time(),"%Y-%m-%d_%H%M%S"),".pdf"),
              height=16, useDingbats=F)
  }
  #
  #-------------------------------------------------------------------------

  {
    # Sort antigens in alphabetical order
    SBA$MFI <- SBA$MFI[, order(colnames(SBA$MFI))]
    SBA$COUNT <- SBA$COUNT[, order(colnames(SBA$COUNT))]
    SBA$BEADS <- SBA$BEADS[order(rownames(SBA$BEADS)),]
    SBA$CT <- SBA$CT[, order(colnames(SBA$CT))]

    # Data transformation
    SBA <- ap_norm2(SBA, MADlimits=seq(0,70,5))
  }

  {
    # plots
    usedata <- SBA
    ap_negbeads(usedata,
                filename=paste0("Results/Second shipment/QC/",analysis,"NegBeads_",format(Sys.time(),"%Y-%m-%d_%H%M%S"),".pdf"),
                width=15, height=10, useDingbats=F)

    ap_agresults(x=usedata,
                 filename=paste0("Results/Second shipment/",analysis,"AntigenResults_",
                                 min(usedata$COKEY$xmad, na.rm=T), "-",max(usedata$COKEY$xmad, na.rm=T),"xMAD_",
                                 format(Sys.time(),"%Y-%m-%d_%H%M%S"),".pdf"),
                 height=15,
                 useDingbats=F)


    # Save data as an excel-file.
    WriteXLS(usedata[-which(names(usedata) %in% c("CT","FILTERINFO","BINARY","DENS","AGCO_CONT",
                                                  "REACTSUM_AG", "REACTFREQ_AG", "REACTSUM_SAMP", "REACTFREQ_SAMP"))],
             paste0("Results/Second shipment/", analysis, "Data_", format(Sys.time(),"%Y-%m-%d_%H%M%S"),".xlsx"),
             row.names=T)
  }
  sessionInfo() # See all version numbers used in the current script run.
} # END All



