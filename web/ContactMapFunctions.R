# Some constants and functions for manipulating and visualising contact maps, both real and predicted
# Shaun Kandathil, UCL

# Requirements
require(bio3d, quietly = T)
require(colorRamps, quietly = T)

# Globals

# distance cutoff in Angstrom for defining a contact
DCUT <- 8.0

# Functions

# not strictly contact-map-related. ASSUMES single-chain PDB.

# Derp there is bio3d::pdbseq(). That selects CA atoms and gets their resids...

# get_seq_from_pdb <- function(pdb, string = TRUE){
#     # if a repeated aas occur, their resno entries should differentiate them. Famous last words?
#     seq_tbl <- unique(pdb$atom[,c('resid','resno')])
#     seq <- aa321(seq_tbl$resid)
#
#     if (string){
#         paste0(seq, collapse = '')
#     } else {
#         seq
#     }
# }

# TODO:
# A function to read in a PFRMAT RR contact table and automatically handle present/absent PFRMAT RR line
# A function to write a con file that works with CMView
# New colour scheme for image_contact_mtx; matlab.like2 doesn't play well with projectors
#   For inter + intra this works well:
#   image_contact_mtx(inter+intra, zlim = c(0,1.5), col = c('black', 'red', 'grey20', 'white'))

select_con_atoms_from_pdb <- function (pdb) {
    # pdb should be an object returned by bio3d::read.pdb
    # This function returns an object of the same type, which has records only for atoms relevant in contact assessment

    # get residue numbers for book-keeping.
    # ASSUMES that there is a single chain
    res_nums <- unique(pdb$atom$resno)

    # select CB for the residues that have them, and CA for Gly
    cb_inds <- atom.select(pdb, elety = "CB")
    gly_ca_inds <-
        atom.select(pdb,
            resid = 'GLY',
            elety = 'CA',
            operator = 'AND')
    con_atoms <-
        combine.select(cb_inds,
            gly_ca_inds,
            operator = 'OR',
            verbose = F)

    # con_atoms should now be an object of type 'select' which defines the subset of atoms we want.
    # of course, we should (at least) check that:
    # 1. there is exactly one atom selected per residue (else we may have AltLocs or other problems)
    # 2. there are the correct number of atoms (== number of residues).
    # Both can be checked somewhat by comparing the residue numbers from the 'trimmed' pdb with res_nums above.

    trim_pdb <- trim.pdb(pdb = pdb, inds = con_atoms)
    # stopifnot(all(trim_pdb$atom$resno == res_nums)) # For a properly concatenated file this should be okay...

    return(trim_pdb)
}

get_con_mtx_from_pdb <- function(pdb, return.logical = TRUE, include.diag = FALSE) {
    # pdb should be an object returned by bio3d::read.pdb

    con_atoms_pdb <- select_con_atoms_from_pdb(pdb)

    xyz <- con_atoms_pdb$atom[,c('x','y','z')]

    d <- as.matrix(dist(xyz))
    dimnames(d) <- NULL

    con <- apply(d, c(1,2), function(x) {x < DCUT}) # contactbench uses < (rather than <= )

    if (!include.diag)
        diag(con) <- FALSE

    if (!return.logical) {
        nres <- dim(d)[1]
        con <- matrix(as.integer(con), nrow = nres)
    }

    return(con)
}

get_dist_mtx_from_pdb <- function(pdb, normalise = F) {
    # pdb should be an object returned by bio3d::read.pdb

    con_atoms_pdb <- select_con_atoms_from_pdb(pdb)

    xyz <- con_atoms_pdb$atom[,c('x','y','z')]

    d <- as.matrix(dist(xyz))

    if(normalise) {
        return(d/max(d))
    } else {
        return (d)
    }
}

get_con_tbl_from_pdb <- function(pdb) {
    mtx <- get_con_mtx_from_pdb(pdb)
    inds <- which(mtx, arr.ind = T)
    rownames(inds) <- NULL
    n <- nrow(inds)
    res <- cbind(inds, rep(0, n), rep(8, n), rep(1, n))
    colnames(res) <- paste0('V', 1:5)
    return(as.data.frame(res))
}

get_dist_tbl_from_pdb <- function(pdb, normalise = T) {
    stop('This function is not ready to use yet.')
    mtx <- get_dist_mtx_from_pdb(pdb, normalise = normalise)
    ## NOT READY YET!
    inds <- which(mtx, arr.ind = T)
    rownames(inds) <- NULL
    n <- nrow(inds)
    res <- cbind(inds, rep(0, n), rep(8, n), rep(1, n)) # <(") (V)_(",,")_(V)
    colnames(res) <- paste0('V', 1:5)
    return(as.data.frame(res))
}

# Convert a 'PFRMAT RR' contact table to a (symmetric) matrix with the contact scores.
pfrmat_to_mtx <-
    function(con_table, nresidues, min_contact_score = 0.0) {
        # first get protein length.
        # By default, ASSUMES that con_table is a full contact list (has data for every possible residue pair, and only the upper/lower triangle)

        # Optional filtering by contact scores here
        if (min_contact_score > 0.0)
            con_table <-
                con_table[which(con_table[, 5] >= min_contact_score), ]

        ncontacts <- nrow(con_table)

        if (missing(nresidues))
            nresidues <- max(con_table[, 1:2])

        mtx <- matrix(0.0, nrow = nresidues, ncol = nresidues)

        # The contacts could be listed in arbitrary order, so I think this has to be the way it is
        # maybe an 'apply' version exists.
        for (i in 1:ncontacts) {
            mtx[con_table[i, 1], con_table[i, 2]] <- con_table[i, 5]
        }

        # These loops could be optimised with R-isms
        for (i in 1:(nresidues - 1)) {
            for (j in (i + 1):nresidues) {
                mtx[j, i] <- mtx[i, j]
            }
        }

        mtx
    }

mtx_to_pfrmat <- function(mtx, min_contact_score = 0.0) {
    #inverse function of pfrmat_to_mtx, with some caveats.
    mtx <- as.matrix(mtx)
    dims <- dim(mtx)
    stopifnot (dims[1] == dims[2])
    nres <- dims[1]

    if (is.logical(mtx))
        mtx <- matrix(as.double(mtx), nrow = nres)

    # need only deal with one triangle of mtx, but make sure it has some valid entries (based on min_contact_score).
    # if not, check the other triangle. If both are empty, issue warning.
    tri <- upper.tri(mtx)

    tri_inds <- which(tri, arr.ind = T)

    con_scores <- mtx[tri_inds]

    # These are row indices for tri_inds, argh
    inds <- which(con_scores > min_contact_score) # behaviour different from other functions that have min_contact_score

    n_con <- length(inds) # This check is not good enough if min_contact_score == 0.0

    # if (n_con == 0) {
    #     tri <- lower.tri(mtx)
    #     tri_inds <- which(tri, arr.ind = T)
    #     con_scores <- mtx[tri_inds]
    #     inds <- which(con_scores > min_contact_score) # behaviour different from other functions that have min_contact_score
    #     n_con <- length(inds)
        if(n_con == 0)
            warning(paste('mtx_to_pfrmat: Zero contacts remain after conversion. min_contact_score =', min_contact_score))
    # }

    #finally get the contact scores
    filtered_tri_inds <- tri_inds[inds,]
    # order the indices for comparability to deepcov output; this really is inconsequential
    # filtered_tri_inds <- filtered_tri_inds[order(filtered_tri_inds[,1], filtered_tri_inds[,2]),]

    selected_contact_scores <- mtx[filtered_tri_inds]

    res <- cbind(filtered_tri_inds, rep(0, n_con), rep(8, n_con), selected_contact_scores)
    colnames(res) <- paste0('V', 1:5)
    return(as.data.frame(res))


}

# aliases for the functions; could be problematic?
mtx_to_tbl <- mtx_to_pfrmat
tbl_to_mtx <- pfrmat_to_mtx

image_contact_mtx <- function(mtx, n_col_levels = 100, col = colorRamps::matlab.like2(n_col_levels), zlim = range(mtx), xlab, ylab, main, ...) {
    # mtx should be a square symmetric matrix returned by pfrmat_to_mtx or similar
    # get nresidues

    d <- dim(mtx)
    stopifnot(length(d) == 2)
    stopifnot(d[1] == d[2])

    nresidues <- d[1]
    ax.lab <- paste('Residue index', c('i', 'j'))
    image(
        x = 1:nresidues,
        y = 1:nresidues,
        z = mtx,
        zlim = zlim,
        col = col,
        xlab = ifelse(missing(xlab), ax.lab[1], xlab),
        ylab = ifelse(missing(ylab), ax.lab[2], ylab),
        main = ifelse(missing(main), 'Contact map', main),
        useRaster = T,
        ...
    )
}

# Get the intersection of two PFRMAT RR-type contact maps.
# optionally apply a contact score cutoff to both maps before calculating the intersection
# con_scores_from is a string, one of 'neither', 'tbl1', 'tbl2', or 'both'.
#   when set to 'neither', all the common contacts are listed with their score set to 1, and the individual contact scores are discarded.
#   when set to either 'tbl1' or 'tbl2', the scores from the specified table are returned.
#   when set to 'both', both sets of contact scores are returned in the table, making the result no longer a strict PFRMAT RR table.
common_contacts_tbl <- function(tbl1, tbl2, min_contact_score = 0.0, con_scores_from='neither') {

    if (!(con_scores_from %in% c('neither', 'tbl1', 'tbl2', 'both')))
        stop("con_scores_from should be one of 'neither', 'tbl1', 'tbl2', or 'both'.")

    # Optional filtering by contact scores here
    if (min_contact_score > 0.0) {
        tbl1 <- tbl1[which(tbl1[, 5] >= min_contact_score), ]
        tbl2 <- tbl2[which(tbl2[, 5] >= min_contact_score), ]
    }

    # merge() with the default of all=FALSE returns an inner join (set intersection) of the two tables
    res <- merge(tbl1, tbl2, by=c(1,2))

    # If res has no entries, it could be that the first two columns are switched in one of the tables
    if (nrow(res) == 0)
        #     tbl1 <- tbl1[,c(2,1,3,4,5)]
        #
        # res <- merge(tbl1, tbl2, by=c(1,2))
        res <- merge(tbl1, tbl2, by.x=c(1,2), by.y=c(2,1))

    # If there are still no entries in res, then output a warning.
    if (nrow(res) == 0)
        warning('No common contacts found in supplied tables at specified cutoff!', immediate. = T)

    if (con_scores_from == 'neither'){
        res <- cbind(res[,c(1,2,3,4)], rep(1, nrow(res)))
        names(res) <- paste0('V', 1:5)
        return(res)
    }

    remove_cols <- c(6,7) # these are always removed

    if (con_scores_from == 'tbl1') {
        remove_cols <- c(remove_cols, 8)
    } else {
        remove_cols <- c(4, remove_cols)
    }

    return(res[, -remove_cols])
}

# Get the intersection of two matrix-type contact maps.
# optionally apply a contact score cutoff to both maps before calculating the intersection
common_contacts_mtx <- function(mtx1, mtx2, min_contact_score = 0.0) {
    stop('This function is not implemented.')
    # Notes: could simply turn into a binary (logical) matrix with something like
    # which((mtx1 > min_contact_score) & (mtx2 > min_contact_score)), arr.ind=T)
    # but it would be nice to have both score sets; that would need a list-like structure
}

filter_short_range_contacts_tbl <- function(tbl, max_co = 4) {
    # max_co = 4 means that ALL contacts with order i --> i+4 and below will be removed

    keep <- abs(tbl[ ,2] - tbl[ ,1]) > max_co
    res <- tbl[keep, ]

    if (nrow(res)==0)
        warning(paste('Zero contacts remain after filtering; max_co was', max_co), immediate. = T, call. = T)

    return(res)
}

plot_common_contacts_tbl <- function(tbl1, tbl2, nresidues = max(tbl1[,1:2], tbl2[,1:2]), min_contact_score = 0.0, con_scores_from='neither') {
    # the default way of setting nresidues ASSUMES that at least one of tbl1 or tbl2 lists every possible contact.

    common <- common_contacts_tbl(tbl1, tbl2, min_contact_score = min_contact_score, con_scores_from = con_scores_from)
    image_contact_mtx(pfrmat_to_mtx(common, nresidues = nresidues))
}

compare_contact_list_to_pdb_tbl <- function(tbl, pdb, min_contact_score = 0.0, min_co = 5, verbose = TRUE) {
    # By default, we only evaluate i --> i+5 or higher order contacts.
    # For only long-range contacts, the convention is min_co = 24.
    # pdb should be an object returned by bio3d::read.pdb
    true_con <- get_con_tbl_from_pdb(pdb)

    if(min_co > 0){
        true_con <- filter_short_range_contacts_tbl(true_con, max_co = (min_co-1))
        tbl <- filter_short_range_contacts_tbl(tbl, max_co = (min_co-1))
    }
    common_con <- common_contacts_tbl(true_con, tbl, con_scores_from = 'tbl2', min_contact_score = min_contact_score)

    ntrue <- nrow(true_con)
    npred <- nrow(tbl)
    ncommon <- nrow(common_con)

    if (verbose){
        writeLines(paste("Percentage of predicted contacts in tbl that are correct (Precision): ", ncommon / npred * 100))
        writeLines(paste("Percentage of true contacts that are in tbl (Coverage): ", ncommon / ntrue * 100))
    }

    return (data.frame(precision = ncommon / npred, coverage = ncommon / ntrue, n_predicted = npred))
}


plot_con <- function(pred_con, min_contact_score=0.5, add=F, ...) {
    if (add)
        pltfn <- points
    else
        pltfn <- plot

    pred_con <- pred_con[which(pred_con$score >= min_contact_score),]
    pltfn(x = pred_con$from, y= pred_con$to, ...)
}


#### EXAMPLE USAGE
# pdb <- read.pdb(file = '11AS_renumbedred_homodimer.pdb')
# pdb_con <- get_con_tbl_from_pdb(pdb = pdb)
# con <- read.table(file = '~/projects/11ASA.con')
#
# ## optional filtering to leave only long-range contacts
# pdb_con <- filter_short_range_contacts_tbl(tbl = pdb_con, max_co = 23)
# con <- filter_short_range_contacts_tbl(tbl = con, max_co = 23)
#
# TPs <- common_contacts_tbl(tbl1 = pdb_con, tbl2 = con, min_contact_score = 0.5, con_scores_from = 'tbl2')
# image_contact_mtx(pfrmat_to_mtx(TPs))

#### EXAMPLE USAGE from .con to .pdf
if (interactive()) {
    con_fname <- read.table(file = '~/Downloads/3d93ae6c-5583-11e9-a887-00163e100d53.mp_stage1')
    pdf_fname <- '~/test.pdf'
} else {
    argv <- commandArgs(trailingOnly = T)
    argc <- length(argv)
    if (argc < 2) {
            stop('ContactMapFunctions.R: usage: Rscript ContactMapFunctions.R <contact file> <PDF filename for output>')
    }
    con_fname <- argv[1]
    pdf_fname <- argv[2]
}

con <- read.table(con_fname)
pdf(file=pdf_fname, height=10, width=10)
image_contact_mtx(pfrmat_to_mtx(con), zlim = c(0,1))
dev.off()

