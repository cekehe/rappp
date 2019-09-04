#### EXAMPLE SCRIPT USING rappp_0.1.0.9004
{ # Run from here for full script
  { # Run from here for all loading, QC, transformation and result files.
    { # Run from here to load and format all input, incl packages

      ################### LOAD Packages ############################################
      {
        # ------------------------ Packages & functions ---------------------------------
        library(rappp) # all ap_-functions
        library(SBA) # read_FlexMAP3D()
        #-----------------------------------------------------------------------
      }

      ################### LOAD DATA & INFO #########################################
      { # Load and format all input
        # ------------------------ Basic Project input ---------------------------------
        {
          analysis <- "DRNA_" # Used as prefix for all filenames, not needed for rappp-functions though.
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
        #
        #-----------------------------------------------------------------------

        # -------- Adjusting matrices ------------------------------------------
        # Fix with gene-names in beads
        # PrESTs mapping to more than one gene are shortened to only the first gene followed by an asterix (*).
        # If the first gene has the structure ACxxxxxx, the second gene name will be used instead.
        # The new gene-names are combined with HPRR to use as an identifier for the rest of the script.
        # You may need to change the spearator from "," to ";", and the column for Gene-names in strsplit(),
        # Check your bead-info.
        {
          beads <- data.frame(beads_org,
                              GeneShort=unlist(lapply(strsplit(as.character(beads_org$Gene.name), ","),
                                                      function(y) ifelse(length(y) > 1 & !grepl("^AC[0-9]{4,}", y[1]),
                                                                         paste0(y[1],"*"),
                                                                         ifelse(length(y) > 1 & grepl("^AC[0-9]{4,}", y[1]),
                                                                                paste0(y[2],"*"), y[1])))),
                              GeneShortAC=unlist(lapply(strsplit(as.character(beads_org$Gene.name), ","),
                                                        function(y) ifelse(length(y) > 1, paste0(y[1],"*"), y[1]))),
                              stringsAsFactors=F)

          beads <- data.frame(Gene_HPRR=ifelse(beads$Type == "PrEST",
                                               paste0(beads$GeneShort,"_HPRR",beads$PrEST.ID),
                                               paste(beads$Antigen.name)), # Column which also includes names for the control beads
                              beads, stringsAsFactors=F)

          rownames(beads) <- beads$Gene_HPRR
        }


        # Rename row- and colnames
        {
          rownames(sampleinfo) <- ifelse(grepl("EMPTY|MIX",sampleinfo$sample_name), # LIMS-IDs
                                         as.character(sampleinfo$tube_label), # Non LIMS-IDs, and with replicate and buffer wells
                                                                              # with conescutive numbering to allow as rownames.
                                         as.character(sampleinfo$sample_name))

          colnames(mfi) <- beads$Gene_HPRR ; rownames(mfi) <- rownames(sampleinfo)
          colnames(count) <- beads$Gene_HPRR ; rownames(count) <- rownames(sampleinfo)


          colnames(mfi_ct) <- beads$Gene_HPRR
        }
        #
        #-----------------------------------------------------------------------

        #--------------- Create data list --------------------------------------------
        #
        # The below structure is used by all ap_-functions used in this script,
        # and is built upon if function output is stored in a new/the same variable
        # (except from ap_reactsummary2 or ap_agresults).
        SBA <- list(MFI=mfi,
                    COUNT=count,
                    SAMPLES=sampleinfo,
                    BEADS=beads,
                    CT=mfi_ct,
                    FILTERINFO=NULL)
        #
        #-------------------------------------------------------------------------
      } # END Load and format all input
    } # END Load and format all input, incl packages


    ################### QC AND DATA TRANSFORMATION #########################################
    { # Run for all QC filtering
    # -------- QC plots and filtering ------------------------------------------
      #
      shouldplot=T
      ## Coupling efficiency test
      SBA <- ap_ct(SBA, grep("Bare", SBA$BEADS$Gene_HPRR), empty_co_multiple=3, shouldplot=shouldplot,
                   filename=paste0("Results/QC/",analysis,"Coupling-efficiency_",format(Sys.time(),"%Y-%m-%d_%H%M%S"),".pdf"),
                   width=30, height=6, useDingbats=F)

      ## Anti-human IgG
      SBA <- ap_igx(SBA, grep("IgG", SBA$BEADS$Gene_HPRR), IgType="G", IgX_cutoff=5000, cosfac=c(3, -3),
                    shouldplot=shouldplot,
                    filename=paste0("Results/QC/",analysis,"ahIgG_",format(Sys.time(),"%Y-%m-%d_%H%M%S"),".pdf"),
                    width=10, height=6, useDingbats=F)

      ## Bead count
      SBA <- ap_count(SBA, labels="Gene_HPRR", protein="GeneShort", agID="PrEST.ID",
                      samp_co=32, bead_flag=32, bead_filter=16, N_filter=0, shouldplot=shouldplot,
                      filename=paste0("Results/QC/",analysis,"BeadCount_",format(Sys.time(),"%Y-%m-%d_%H%M%S"),".pdf"),
                      width=12, height=10, useDingbats=F)

      ## Manual annotations
      SBA$SAMPLES$Filtered[grep("rep|pool|commercial|mix|blank|empty|buffer", SBA$SAMPLES$sample_name, ignore.case=T)] <- "Control"
      SBA$BEADS$Filtered[grep("ebna|ahig|human", SBA$BEADS$Gene_HPRR, ignore.case=T)] <- "PosControl"
      SBA$BEADS$Filtered[grep("bare|empty|his6abp|hisabp", SBA$BEADS$Gene_HPRR, ignore.case=T)] <- "NegControl"
      #
      #-------------------------------------------------------------------------
    } # END QC filtering

    { # Run for all data transformation
      # -------- Data transformation ------------------------------------------
      #
      # Sort antigens in alphabetical order
      SBA$MFI <- SBA$MFI[, order(colnames(SBA$MFI))]
      SBA$COUNT <- SBA$COUNT[, order(colnames(SBA$COUNT))]
      SBA$BEADS <- SBA$BEADS[order(rownames(SBA$BEADS)), ]
      SBA$CT <- SBA$CT[, order(colnames(SBA$CT))]

      # Data transformation
      SBA <- ap_norm2(SBA,
                      MADlimits = seq(0,70,5),
                      na.rm = TRUE,
                      check.names = FALSE,
                      mad_constant = 1,
                      mad_low = FALSE,
                      mad_high = FALSE,
                      score_rightmost.closed = FALSE,
                      score_left.open = FALSE,
                      score_all.inside = FALSE,
                      coselect_slope_cutoff = -0.5,
                      coselect_offset = 0.1,
                      coselect_bw = 0.1)
      #
      #-------------------------------------------------------------------------
    } # END data transformation

    ################### RESULTS #########################################
    { # Run for results and data output
      # -------- Reactivity results and data output ----------------------------
      #
      ## Calculate reactivity frequencies per sample group and plot it
      react_summary <- ap_agresults(x=SBA,
                                    samplegroups = NULL,
                                    filename=paste0("Results/",analysis,"AntigenResults_",
                                                    min(SBA$CUTOFF_KEY$xmad, na.rm=T), "-",
                                                    max(SBA$CUTOFF_KEY$xmad, na.rm=T),"xMAD_",
                                                    format(Sys.time(),"%Y-%m-%d_%H%M%S"),".pdf"),
                                    height=15,
                                    useDingbats=F)

      ## Calculate reactivity frequencies per sample group without plotting it (same data output as if plotting)
      react_summary <- ap_reactsummary2(SBA) # Without plot

      ## Save data as an excel-file.
      ap_excel(append(SBA, react_summary),
               elements = c("MFI", "MADS", "SCORE", "BINARY_CO",
                            "REACTIVITIES", "ANTIGEN_CUTOFFS", "CUTOFF_KEY",
                            "ANTIGENS", "SAMPLES", "COUNT"),
               filename = paste0("Results/", analysis, "Data_",
                                 format(Sys.time(),"%Y-%m-%d_%H%M%S"),".xlsx"),
               row.names = TRUE)

      ## Print out summaries regarding number of Uniprots, fragment per uniprot, aa-lengths etc.
      ap_summary(SBA)

      ## Image for cutoff key
      ap_cutoffs2image(cutoffkey = SBA$CUTOFF_KEY)
      #
      #-------------------------------------------------------------------------
    } # END Results and data output
  } # END all loading, QC, transformation and results output.

  # -------- Other plots to get an overview of the data --------------------
  #
  {
    ## Signal overview
    ap_overview(SBA, filename=paste0("Results/QC/",analysis,"Signal-overview_",format(Sys.time(),"%Y-%m-%d_%H%M%S"),".pdf"),
                width=20, height=15, useDingbats=F, cex.axis=0.3)

    ## CV and replicates
    ap_rep(SBA, iter=500, filename=paste0("Results/QC/",analysis,"reproducibility_",format(Sys.time(),"%Y-%m-%d_%H%M%S"),".pdf"),
           width=12, height=12, useDingbats=F)

    ## tSNE
    # Note that this function doesn't use the the whole SBA-list as input,
    # you only have the data you want to plot as input.
    tsne_perp(SBA$MFI, perp=c(2,5,10,50), sqrt=T, iterations=1000,
              groups=as.factor(SBA$SAMPLES$AssayNum),
              names=rownames(SBA$MFI),
              legend=T,
              filename=paste0("Results/QC/",analysis,"tSNE_",format(Sys.time(),"%Y-%m-%d_%H%M%S"),".pdf"),
              height=16, useDingbats=F)

    ## Scoring and the neg-control beads
    ap_negbeads(SBA,
                filename=paste0("Results/QC/",analysis,"NegBeads_",format(Sys.time(),"%Y-%m-%d_%H%M%S"),".pdf"),
                width=15, height=10, useDingbats=F)
  }
  #
  #-------------------------------------------------------------------------

  sessionInfo() # See all version numbers used in the current script run.
} # END All
