/* Copyright (c) 2023, Huiyu Li and Ao Ma, University of Illinois at Chicago.
 * This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
* This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
 * You should have received a copy of the GNU General Public License along with this program. If not, see <https://www.gnu.org/licenses/>.
 */

/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2013,2014,2015,2016,2017,2018, by the GROMACS development team, led by
 * Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
 * and including many others, as listed in the AUTHORS file in the
 * top-level source directory and at http://www.gromacs.org.
 *
 * GROMACS is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 *
 * GROMACS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with GROMACS; if not, see
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at http://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org.
 */
#include "gmxpre.h"

#include "mdoutf.h"

#include "gromacs/commandline/filenm.h"
#include "gromacs/domdec/collect.h"
#include "gromacs/domdec/domdec_struct.h"
#include "gromacs/fileio/checkpoint.h"
#include "gromacs/fileio/gmxfio.h"
#include "gromacs/fileio/tngio.h"
#include "gromacs/fileio/trrio.h"
#include "gromacs/fileio/xtcio.h"
#include "gromacs/fileio/xvgr.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdlib/mdrun.h"
#include "gromacs/mdlib/trajectory_writing.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/imdoutputprovider.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/mdtypes/state.h"
#include "gromacs/timing/wallcycle.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/pleasecite.h"
#include "gromacs/utility/smalloc.h"

/* added by Huiyu Li */
#include "gromacs/mdtypes/inputrec.h"
/* end of adding by Huiyu Li*/

struct gmx_mdoutf {
    t_fileio               *fp_trn;
    t_fileio               *fp_xtc;
    gmx_tng_trajectory_t    tng;
    gmx_tng_trajectory_t    tng_low_prec;
    int                     x_compression_precision; /* only used by XTC output */
    ener_file_t             fp_ene;
    const char             *fn_cpt;
    gmx_bool                bKeepAndNumCPT;
    int                     eIntegrator;
    gmx_bool                bExpanded;
    int                     elamstats;
    int                     simulation_part;
    FILE                   *fp_dhdl;
    int                     natoms_global;
    int                     natoms_x_compressed;
    gmx_groups_t           *groups; /* for compressed position writing */
    gmx_wallcycle_t         wcycle;
    rvec                   *f_global;
    gmx::IMDOutputProvider *outputProvider;


    /* added by Huiyu Li, output protein atoms of system into a new trr file */
    gmx_bool               selftrr;   /* 1 for output self trr, 0 for not */
    t_fileio              *fpst;      /* file pointer for self trr file */
    int                    selfnatom;
    int                   *hlstep;
    int                    deffreq;  /* nstxout, nstvout, nstfout in mdp file */
                                     /*  will control new trr output frequency */
                                     /* deffreq will control Gromacs trr output frequency */
    /* end of adding by Huiyu Li */
};

/* added by Huiyu Li, output self trr file */
void init_self_mdoutf(gmx_mdoutf_t *outf, t_inputrec *ir, int *hlstep, char *dir)
{
  char nm[800];

  fprintf(stderr, "\nHL initialize self trr file\n\n");
  if (ir->bias.selftrr == 1){
    (*outf)->selftrr = TRUE;
  } else if (ir->bias.selftrr == 0){
    (*outf)->selftrr = FALSE;
  } else{
    fprintf(stderr, "\nHL Error: init_self_mdoutf wrong ir->bias.selftrr %5d\n\n", ir->bias.selftrr);
    exit(1);
  }
  (*outf)->selfnatom = ir->bias.selfnatom;
  (*outf)->deffreq = ir->bias.deffreq;
  (*outf)->hlstep = hlstep;

  if ((*outf)->selftrr == 1){
    snprintf(nm, 800, "%s/protein.trr", dir);
    (*outf)->fpst = gmx_trr_open(nm, "w");
  }  
}

char* hl_parse_currentdir(gmx_mdoutf_t *outf)
{
  int i;
  char* curdir;
  char* nulldir;
  size_t size;
  const char *loc = (*outf)->fn_cpt;

  size = strlen(loc);

  for (i = (size - 1); (i >= 0) && (loc[i] != '/'); --i);

  if (loc[i] == '/'){
    curdir = (char*) malloc((i + 1) * sizeof(char));
    if (curdir != NULL){
      strncpy(curdir, loc, (size_t) i);
      curdir[i] = '\0';
      return curdir;
    } else{
      fprintf(stderr,"Memory allocation error.");
    }
  } else{
    fprintf(stderr,"There was no DIR_SEPARATOR in LOC.");
  }

  nulldir = (char*) malloc(800 * sizeof(char));
  snprintf(nulldir, 800, "./");
  return  nulldir;
}

/* end of adding by Huiyu Li */

gmx_mdoutf_t init_mdoutf(FILE *fplog, int nfile, const t_filenm fnm[],
                         const MdrunOptions &mdrunOptions,
                         const t_commrec *cr,
                         gmx::IMDOutputProvider *outputProvider,
                         const t_inputrec *ir, gmx_mtop_t *top_global,
                         const gmx_output_env_t *oenv, gmx_wallcycle_t wcycle)
{
    gmx_mdoutf_t   of;
    const char    *appendMode = "a+", *writeMode = "w+", *filemode;
    gmx_bool       bAppendFiles, bCiteTng = FALSE;
    int            i;

    snew(of, 1);

    of->fp_trn       = nullptr;
    of->fp_ene       = nullptr;
    of->fp_xtc       = nullptr;
    of->tng          = nullptr;
    of->tng_low_prec = nullptr;
    of->fp_dhdl      = nullptr;

    of->eIntegrator             = ir->eI;
    of->bExpanded               = ir->bExpanded;
    of->elamstats               = ir->expandedvals->elamstats;
    of->simulation_part         = ir->simulation_part;
    of->x_compression_precision = static_cast<int>(ir->x_compression_precision);
    of->wcycle                  = wcycle;
    of->f_global                = nullptr;
    of->outputProvider          = outputProvider;

    if (MASTER(cr))
    {
        bAppendFiles = mdrunOptions.continuationOptions.appendFiles;

        of->bKeepAndNumCPT = mdrunOptions.checkpointOptions.keepAndNumberCheckpointFiles;

        filemode = bAppendFiles ? appendMode : writeMode;

        if (EI_DYNAMICS(ir->eI) &&
            ir->nstxout_compressed > 0)
        {
            const char *filename;
            filename = ftp2fn(efCOMPRESSED, nfile, fnm);
            switch (fn2ftp(filename))
            {
                case efXTC:
                    of->fp_xtc                  = open_xtc(filename, filemode);
                    break;
                case efTNG:
                    gmx_tng_open(filename, filemode[0], &of->tng_low_prec);
                    if (filemode[0] == 'w')
                    {
                        gmx_tng_prepare_low_prec_writing(of->tng_low_prec, top_global, ir);
                    }
                    bCiteTng = TRUE;
                    break;
                default:
                    gmx_incons("Invalid reduced precision file format");
            }
        }
        if ((EI_DYNAMICS(ir->eI) || EI_ENERGY_MINIMIZATION(ir->eI)) &&
            (!GMX_FAHCORE &&
             !(EI_DYNAMICS(ir->eI) &&
               ir->nstxout == 0 &&
               ir->nstvout == 0 &&
               ir->nstfout == 0)
            )
            )
        {
            const char *filename;
            filename = ftp2fn(efTRN, nfile, fnm);
            switch (fn2ftp(filename))
            {
                case efTRR:
                case efTRN:
                    /* If there is no uncompressed coordinate output and
                       there is compressed TNG output write forces
                       and/or velocities to the TNG file instead. */
                    if (ir->nstxout != 0 || ir->nstxout_compressed == 0 ||
                        !of->tng_low_prec)
                    {
                        of->fp_trn = gmx_trr_open(filename, filemode);
                    }
                    break;
                case efTNG:
                    gmx_tng_open(filename, filemode[0], &of->tng);
                    if (filemode[0] == 'w')
                    {
                        gmx_tng_prepare_md_writing(of->tng, top_global, ir);
                    }
                    bCiteTng = TRUE;
                    break;
                default:
                    gmx_incons("Invalid full precision file format");
            }
        }
        if (EI_DYNAMICS(ir->eI) || EI_ENERGY_MINIMIZATION(ir->eI))
        {
            of->fp_ene = open_enx(ftp2fn(efEDR, nfile, fnm), filemode);
        }
        of->fn_cpt = opt2fn("-cpo", nfile, fnm);

        if ((ir->efep != efepNO || ir->bSimTemp) && ir->fepvals->nstdhdl > 0 &&
            (ir->fepvals->separate_dhdl_file == esepdhdlfileYES ) &&
            EI_DYNAMICS(ir->eI))
        {
            if (bAppendFiles)
            {
                of->fp_dhdl = gmx_fio_fopen(opt2fn("-dhdl", nfile, fnm), filemode);
            }
            else
            {
                of->fp_dhdl = open_dhdl(opt2fn("-dhdl", nfile, fnm), ir, oenv);
            }
        }

        outputProvider->initOutput(fplog, nfile, fnm, bAppendFiles, oenv);

        /* Set up atom counts so they can be passed to actual
           trajectory-writing routines later. Also, XTC writing needs
           to know what (and how many) atoms might be in the XTC
           groups, and how to look up later which ones they are. */
        of->natoms_global       = top_global->natoms;
        of->groups              = &top_global->groups;
        of->natoms_x_compressed = 0;
        for (i = 0; (i < top_global->natoms); i++)
        {
            if (getGroupType(of->groups, egcCompressedX, i) == 0)
            {
                of->natoms_x_compressed++;
            }
        }

        if (ir->nstfout && DOMAINDECOMP(cr))
        {
            snew(of->f_global, top_global->natoms);
        }
    }

    if (bCiteTng)
    {
        please_cite(fplog, "Lundborg2014");
    }

    return of;
}

ener_file_t mdoutf_get_fp_ene(gmx_mdoutf_t of)
{
    return of->fp_ene;
}

FILE *mdoutf_get_fp_dhdl(gmx_mdoutf_t of)
{
    return of->fp_dhdl;
}

gmx_wallcycle_t mdoutf_get_wcycle(gmx_mdoutf_t of)
{
    return of->wcycle;
}

void mdoutf_write_to_trajectory_files(FILE *fplog, const t_commrec *cr,
                                      gmx_mdoutf_t of,
                                      int mdof_flags,
                                      gmx_mtop_t *top_global,
                                      int64_t step, double t,
                                      t_state *state_local, t_state *state_global,
                                      ObservablesHistory *observablesHistory,
                                      gmx::ArrayRef<gmx::RVec> f_local)
{
    /* added by Huiyu Li */
    gmx_bool whedefout;
    int backupnatom;
    /* end of adding by Huiyu Li */

    rvec *f_global;

    if (DOMAINDECOMP(cr))
    {
        if (mdof_flags & MDOF_CPT)
        {
            dd_collect_state(cr->dd, state_local, state_global);
        }
        else
        {
            if (mdof_flags & (MDOF_X | MDOF_X_COMPRESSED))
            {
                gmx::ArrayRef<gmx::RVec> globalXRef = MASTER(cr) ? makeArrayRef(state_global->x) : gmx::EmptyArrayRef();
                dd_collect_vec(cr->dd, state_local, makeArrayRef(state_local->x), globalXRef);
            }
            if (mdof_flags & MDOF_V)
            {
                gmx::ArrayRef<gmx::RVec> globalVRef = MASTER(cr) ? makeArrayRef(state_global->v) : gmx::EmptyArrayRef();
                dd_collect_vec(cr->dd, state_local, makeArrayRef(state_local->v), globalVRef);
            }
        }
        f_global = of->f_global;
        if (mdof_flags & MDOF_F)
        {
            dd_collect_vec(cr->dd, state_local, f_local, gmx::arrayRefFromArray(reinterpret_cast<gmx::RVec *>(f_global), f_local.size()));
        }
    }
    else
    {
        /* We have the whole state locally: copy the local state pointer */
        state_global = state_local;

        f_global     = as_rvec_array(f_local.data());
    }

    if (MASTER(cr))
    {
        if (mdof_flags & MDOF_CPT)
        {
            fflush_tng(of->tng);
            fflush_tng(of->tng_low_prec);
            ivec one_ivec = { 1, 1, 1 };
            write_checkpoint(of->fn_cpt, of->bKeepAndNumCPT,
                             fplog, cr,
                             DOMAINDECOMP(cr) ? cr->dd->nc : one_ivec,
                             DOMAINDECOMP(cr) ? cr->dd->nnodes : cr->nnodes,
                             of->eIntegrator, of->simulation_part,
                             of->bExpanded, of->elamstats, step, t,
                             state_global, observablesHistory);
        }

        if (mdof_flags & (MDOF_X | MDOF_V | MDOF_F))
        {
            const rvec *x = (mdof_flags & MDOF_X) ? state_global->x.rvec_array() : nullptr;
            const rvec *v = (mdof_flags & MDOF_V) ? state_global->v.rvec_array() : nullptr;
            const rvec *f = (mdof_flags & MDOF_F) ? f_global : nullptr;

            if (of->fp_trn)
            {
              /* added by Huiyu Li, output control for Gromacs trr */
              if (of->selftrr){
                if ((of->hlstep[0] % of->deffreq) == 0){
                  whedefout = true;
                } else{
                  whedefout = false;
                }
              } else{
                whedefout = true;
              }
              /* end of adding by Huiyu Li */

              if (whedefout){  /* this line is added by Huiyu Li, control Gromacs trr output */
                gmx_trr_write_frame(of->fp_trn, step, t, state_local->lambda[efptFEP],
                                    state_local->box, top_global->natoms,
                                    x, v, f);
                if (gmx_fio_flush(of->fp_trn) != 0)
                {
                    gmx_file("Cannot write trajectory; maybe you are out of disk space?");
                }
              } /* this line added by Huiyu Li, Gromacs trr output*/
            }

            /* added by Huiyu Li, output selftrr */
            if (of->selftrr){
              backupnatom = top_global->natoms;
              top_global->natoms = of->selfnatom;
              gmx_trr_write_frame(of->fpst, step, t, state_local->lambda[efptFEP],
                                  state_local->box, top_global->natoms,
                                  x, v, f);
              top_global->natoms = backupnatom;
            }
            /* end of adding by Huiyu Li */


            /* If a TNG file is open for uncompressed coordinate output also write
               velocities and forces to it. */
            else if (of->tng)
            {
                gmx_fwrite_tng(of->tng, FALSE, step, t, state_local->lambda[efptFEP],
                               state_local->box,
                               top_global->natoms,
                               x, v, f);
            }
            /* If only a TNG file is open for compressed coordinate output (no uncompressed
               coordinate output) also write forces and velocities to it. */
            else if (of->tng_low_prec)
            {
                gmx_fwrite_tng(of->tng_low_prec, FALSE, step, t, state_local->lambda[efptFEP],
                               state_local->box,
                               top_global->natoms,
                               x, v, f);
            }
        }
        if (mdof_flags & MDOF_X_COMPRESSED)
        {
            rvec *xxtc = nullptr;

            if (of->natoms_x_compressed == of->natoms_global)
            {
                /* We are writing the positions of all of the atoms to
                   the compressed output */
                xxtc = state_global->x.rvec_array();
            }
            else
            {
                /* We are writing the positions of only a subset of
                   the atoms to the compressed output, so we have to
                   make a copy of the subset of coordinates. */
                int i, j;

                snew(xxtc, of->natoms_x_compressed);
                auto x = makeArrayRef(state_global->x);
                for (i = 0, j = 0; (i < of->natoms_global); i++)
                {
                    if (getGroupType(of->groups, egcCompressedX, i) == 0)
                    {
                        copy_rvec(x[i], xxtc[j++]);
                    }
                }
            }
            if (write_xtc(of->fp_xtc, of->natoms_x_compressed, step, t,
                          state_local->box, xxtc, of->x_compression_precision) == 0)
            {
                gmx_fatal(FARGS,
                          "XTC error. This indicates you are out of disk space, or a "
                          "simulation with major instabilities resulting in coordinates "
                          "that are NaN or too large to be represented in the XTC format.\n");
            }
            gmx_fwrite_tng(of->tng_low_prec,
                           TRUE,
                           step,
                           t,
                           state_local->lambda[efptFEP],
                           state_local->box,
                           of->natoms_x_compressed,
                           xxtc,
                           nullptr,
                           nullptr);
            if (of->natoms_x_compressed != of->natoms_global)
            {
                sfree(xxtc);
            }
        }
        if (mdof_flags & (MDOF_BOX | MDOF_LAMBDA) && !(mdof_flags & (MDOF_X | MDOF_V | MDOF_F)) )
        {
            if (of->tng)
            {
                real  lambda = -1;
                rvec *box    = nullptr;
                if (mdof_flags & MDOF_BOX)
                {
                    box = state_local->box;
                }
                if (mdof_flags & MDOF_LAMBDA)
                {
                    lambda = state_local->lambda[efptFEP];
                }
                gmx_fwrite_tng(of->tng, FALSE, step, t, lambda,
                               box, top_global->natoms,
                               nullptr, nullptr, nullptr);
            }
        }
        if (mdof_flags & (MDOF_BOX_COMPRESSED | MDOF_LAMBDA_COMPRESSED) && !(mdof_flags & (MDOF_X_COMPRESSED)) )
        {
            if (of->tng_low_prec)
            {
                real  lambda = -1;
                rvec *box    = nullptr;
                if (mdof_flags & MDOF_BOX_COMPRESSED)
                {
                    box = state_local->box;
                }
                if (mdof_flags & MDOF_LAMBDA_COMPRESSED)
                {
                    lambda = state_local->lambda[efptFEP];
                }
                gmx_fwrite_tng(of->tng_low_prec, FALSE, step, t, lambda,
                               box, top_global->natoms,
                               nullptr, nullptr, nullptr);
            }
        }
    }
}

void mdoutf_tng_close(gmx_mdoutf_t of)
{
    if (of->tng || of->tng_low_prec)
    {
        wallcycle_start(of->wcycle, ewcTRAJ);
        gmx_tng_close(&of->tng);
        gmx_tng_close(&of->tng_low_prec);
        wallcycle_stop(of->wcycle, ewcTRAJ);
    }
}

void done_mdoutf(gmx_mdoutf_t of)
{
    if (of->fp_ene != nullptr)
    {
        close_enx(of->fp_ene);
    }
    if (of->fp_xtc)
    {
        close_xtc(of->fp_xtc);
    }
    if (of->fp_trn)
    {
        gmx_trr_close(of->fp_trn);
    }
    if (of->fp_dhdl != nullptr)
    {
        gmx_fio_fclose(of->fp_dhdl);
    }
    of->outputProvider->finishOutput();
    if (of->f_global != nullptr)
    {
        sfree(of->f_global);
    }

    gmx_tng_close(&of->tng);
    gmx_tng_close(&of->tng_low_prec);

    sfree(of);


    /* added by Huiyu Li, close self output trr */
    if (of->selftrr){
      fprintf(stderr, "\nHL close self trr file\n");
      gmx_trr_close(of->fpst);
      fprintf(stderr, "HL done closing self trr file\n\n");
    }
}

int mdoutf_get_tng_box_output_interval(gmx_mdoutf_t of)
{
    if (of->tng)
    {
        return gmx_tng_get_box_output_interval(of->tng);
    }
    return 0;
}

int mdoutf_get_tng_lambda_output_interval(gmx_mdoutf_t of)
{
    if (of->tng)
    {
        return gmx_tng_get_lambda_output_interval(of->tng);
    }
    return 0;
}

int mdoutf_get_tng_compressed_box_output_interval(gmx_mdoutf_t of)
{
    if (of->tng_low_prec)
    {
        return gmx_tng_get_box_output_interval(of->tng_low_prec);
    }
    return 0;
}

int mdoutf_get_tng_compressed_lambda_output_interval(gmx_mdoutf_t of)
{
    if (of->tng_low_prec)
    {
        return gmx_tng_get_lambda_output_interval(of->tng_low_prec);
    }
    return 0;
}


/* added by Huiyu Li */

/* set vector elements to zero */
void setzeroHL(rvec *f, int len)
{
  int i, j;
  for (i = 0; i < len; i++){
    for (j = 0; j < 3; j++){
      f[i][j] = 0.0;
    }
  }
}

/* Initialize BIAS struct for collective variable case */
void initializeCV(BIAS *bias)
{
  int i;

  fprintf(stderr, "HL initialize collective variable\n");
  bias->cvcurrbin = -1;
  bias->cvfastbin = -1;

  /* bias bins for collective variable */
  bias->cvbinsize = (bias->cvup - bias->cvlow) / bias->cvnumbin;
  fprintf(stderr, "  HL %5d bins in range [%lf %lf] width %lf\n", 
          bias->cvnumbin, bias->cvlow, bias->cvup, bias->cvbinsize);
  for (i = 0; i <= bias->cvnumbin; i++){
    bias->cvbd[i] = bias->cvlow + i * bias->cvbinsize;
    fprintf(stderr, "    cv bin bd %5d %lf\n", i, bias->cvbd[i]);
  }
  
  /* output coefficient for each dihedral in CV */ 
  for (i = 0; i < bias->numic; i++){
    fprintf(stderr, "  HL CV %5d coeff %lf target %lf\n", i, bias->bic[i].cvcoeff, bias->bic[i].cvcenter);
  }
}

/* malloc and read atom list of BBdih involved in U coordinates */
int **initU(int numdih, char *nm)
{
  int **mat;
  int i, j;
  FILE *fp;
  int one = 4; /* each dihedral has 4 atoms */

  mat = (int **) malloc(numdih * sizeof(int*));
  for (i = 0; i < numdih; i++){
    mat[i] = (int *) malloc(one * sizeof(int));
  }

  fp = fopen(nm, "r");
  for (i = 0; i < numdih; i++){
    for (j = 0; j < one; j++){
      fscanf(fp, "%d", &mat[i][j]);
      mat[i][j] -= 1;   /* atom index need to starts from 0 for C code */
    }
  }
  fclose(fp);

  return mat;
}

/* Initialize BIAS struct, including generating bias bins */
void initializeBIAS(t_inputrec *ir, int natom, char *dir)
{ 
  BIAS *bias;
  char nm[1000];

  fprintf(stderr, "\nHL initialize BIAS struct\n");
  bias = &(ir->bias); 

  /* some quantities */
  bias->everbasin = 0;
  bias->whereach = 0;
  bias->bLastStep = FALSE;
  bias->natom = natom;
  snew(bias->biasf, bias->natom);

  /* open output files */
  snprintf(nm, 1000, "%s/bias_res.dat", dir);
  fprintf(stderr, "HL Output bias information to %s\n", nm);
  bias->fpbias = fopen(nm, "w");

  snprintf(nm, 1000, "%s/bias_force.dat", dir);
  fprintf(stderr, "HL Output bias force to %s\n", nm);
  bias->fpbf = fopen(nm, "w");

  snprintf(nm, 1000, "%s/basin.dat", dir);
  fprintf(stderr, "HL Ouput basin information to %s\n", nm);
  bias->fpbasin = fopen(nm, "w");

  /* initialize bias for different whecv */
  switch (bias->whecv){
    case 1:
      /* collective variable */
      initializeCV(bias);
      break;
    default:
      fprintf(stderr, "\nHL error: unknown whecv %d in initializeBIAS\n\n", bias->whecv);
      exit(1);
  }

  /* calculate BBdih involved in U vector */
  if (bias->wheU > 0){
    snprintf(nm, 1000, "%s/bbdih.dat", dir);
    fprintf(stderr, "HL Output BBdih to %s\n", nm);
    bias->fpdih = fopen(nm, "w");
    bias->ubbdih = (double *) malloc((bias->wheU) * sizeof(double));
    bias->uatom = initU(bias->wheU, bias->nmu);
    fprintf(stderr, "  uatom %5d %5d %5d %5d\n       %5d %5d %5d %5d\n", bias->uatom[0][0], bias->uatom[0][1], bias->uatom[0][2], bias->uatom[0][3], bias->uatom[bias->wheU-1][0], bias->uatom[bias->wheU-1][1], bias->uatom[bias->wheU-1][2], bias->uatom[bias->wheU-1][3]);
  }

  fprintf(stderr, "HL finish initialization of BIAS\n\n");
}

/* close files and free malloc space */
void freeBIAS(t_inputrec *ir)
{
  BIAS *bias;
  int i;
 
  fprintf(stderr, "\nHL freeBIAS\n");
 
  bias = &(ir->bias);
  sfree(bias->biasf);
  fclose(bias->fpbias);
  fclose(bias->fpbf);
  fclose(bias->fpbasin);

  if (bias->wheU > 0){
    fclose(bias->fpdih);
    free(bias->ubbdih);
    for (i = 0; i < bias->wheU; i++){
      free(bias->uatom[i]);
    }
    free(bias->uatom);
  }

  fprintf(stderr, "\nHL free boundaries\n");
  free(bias->cvbd);
  for (i = 0; i < bias->numic; i++){
    free(bias->bic[i].bd);
  }
  sfree(bias->bic);

  fprintf(stderr, "\nHL done freeBIAS\n==========\n\n");
}

/* calculate dihedral */
int rvec_sub_pbcHL(const rvec a, const rvec b, rvec dx, rvec boxsize)
{
  int p, i;
  rvec_sub(a, b, dx);

  p = 0;
  for (i = 0; i < 3; i++){
    if (dx[i] > 0.5 * boxsize[i]){
      dx[i] -= boxsize[i];
      p = 1;
    } else if (dx[i] < -0.5 * boxsize[i]){
      dx[i] += boxsize[i];
      p = -1;
    }
  }

  return p;
}

real dih_angle_HL(const rvec xi,const rvec xj,const rvec xk,const rvec xl, rvec boxsize)
{
  real ipr,phi;
  rvec r_ij, r_kj, r_kl, m, n;
  real sign;
  int p1,p2,p3;

  p1 = rvec_sub_pbcHL(xi, xj, r_ij, boxsize);     /* In bonded.cpp, dih_angle() function, they use pbc_rvec_sub */
  p2 = rvec_sub_pbcHL(xk, xj, r_kj, boxsize);
  p3 = rvec_sub_pbcHL(xk, xl, r_kl, boxsize);

  cprod(r_ij,r_kj,m);                   /*  9           */
  cprod(r_kj,r_kl,n);                   /*  9           */
  phi=gmx_angle(m,n);                   /* 49 (assuming 25 for atan2) */
  ipr=iprod(r_ij,n);                    /*  5           */
  sign=(ipr<0.0)?-1.0:1.0;
  phi=sign*phi;                         /*  1           */
                                        /* 82 TOTAL     */

/* ** This is for testing purpose only *********
  if ((p1!=0)||(p2!=0)||(p3!=0)){
    fprintf(stderr, "\nError: p %5d %5d %5d phi %lf\n\n", p1,p2,p3,phi);
    exit(1);
  } else{
    fprintf(stderr, " correct %lf\n", phi);
  }
 */

  return phi;
}

/* calculate dihedral for one BIC */
void calcdihBIC(const rvec x[], BIC *bic, rvec boxsize)
{
  int ai, aj, ak, al;

  ai = bic->atom[0];
  aj = bic->atom[1];
  ak = bic->atom[2];
  al = bic->atom[3];
  bic->val = dih_angle_HL(x[ai], x[aj], x[ak], x[al], boxsize);
}

/* check in which bin, the val locates */
int whichbiasbin(double val, int numbin, double *bd)
{
  int i;

  for (i = 0; i < numbin; i++){
    if ((val >= bd[i]) && (val < bd[i+1])){
      return i;
    }
  }

  /* return -100, if the value is on the left hand side of the range */
  if (val < bd[0]){
    return -100;
  }

  /* return -1, if the value is on the right hand side of the range */
  if (val >= bd[numbin]){
    return -1;
  }

  fprintf(stderr, "\nHL error: whichbiasbin %lf for %5d bins in [%lf %lf]\n\n",
          val, numbin, bd[0], bd[numbin]);
  exit(1);
}

double calcCV(int numic, BIC *bic, int gbf, int signU)
{
  /* calculate value of collective variable */
  /* value of each dihedral must be calculated before this */
  int i;
  double cv;

  cv = 0.0;
  switch (gbf){
    case 1:
    case 3:
      /* collective variable is sum of cosine functions */
      for (i = 0; i < numic; i++){
        cv += bic[i].cvcoeff * (bic[i].cvshift + cos(bic[i].val - bic[i].cvcenter));
      }
      break;
    case 2:
      /* collective variable is sum of dihedrals */
      for (i = 0; i < numic; i++){
        cv += bic[i].cvcoeff * (bic[i].cvshift + (bic[i].val - bic[i].cvcenter));
      }
      break;
    case 4:
      /* collective variable is sum of cosine function with absolute coefficient */
      for (i = 0; i < numic; i++){
        cv += fabs(bic[i].cvcoeff) * (bic[i].cvshift + cos(bic[i].val - bic[i].cvcenter));
      }
      break;
    case 5:
      /* see case 5~8 on word note on Jul 10 2022*/
      for (i = 0; i < numic; i++){
        cv += signU * bic[i].cvcoeff * bic[i].val;
      }
      break;
    case 6:
      for (i = 0; i < numic; i++){
        cv += signU * fabs(bic[i].cvcoeff) * cos(bic[i].val - bic[i].cvcenter);
      }
      break;
    case 7:
      for (i = 0; i < numic; i++){
        cv += signU * bic[i].cvcoeff * cos(bic[i].val - bic[i].cvcenter);
      }
      break;
    default:
      fprintf(stderr, "\nError: unknown gbf %d in calcCV\n\n", gbf);
      exit(1);
  }

  return cv;
}

/* calculate generalized bias force */
double calc_gene_biasForce(BIAS *bias, double deltaCV, int gbf, int indexbic)
{
  double geneF;

  geneF = 0.0;
  switch (gbf){
    case 1:
      /* fprintf(stderr, "gbf %d CV is sum of cosine functions\n", gbf); */
      geneF = -1.0 * bias->cvstrength * bias->bic[indexbic].cvcoeff * deltaCV *
              sin(bias->bic[indexbic].val - bias->bic[indexbic].cvcenter);
      break;
    case 2:
      /* fprintf(stderr, "gbf %d CV is sum of internal coordinates them selves\n", gbf); */
      geneF = 1.0 * bias->cvstrength * bias->bic[indexbic].cvcoeff * deltaCV;
      break;
    case 3:
      /* fprintf(stderr, "gbf %d CV is sum of cosine, but geneF devided by -sine\n", gbf); */
      geneF = 1.0 * bias->cvstrength * bias->bic[indexbic].cvcoeff * deltaCV;
      break;
    case 4:
      /* case 4: fprintf(stderr, "gbf %d CV is sum of cosine with solute coefficient, but geneF divided by -sine\n", gbf)*/
      geneF = 1.0 * bias->cvstrength * bias->bic[indexbic].cvcoeff * deltaCV;
      break;
    case 5:
      /* see case 5~8 on word note on Jul 10 2022 */
      geneF = 1.0 * bias->cvstrength * deltaCV * bias->signU * bias->bic[indexbic].cvcoeff;
      break;
    case 6:
      geneF = -1.0 * bias->cvstrength * deltaCV * bias->signU * fabs(bias->bic[indexbic].cvcoeff)
              * sin(bias->bic[indexbic].val - bias->bic[indexbic].cvcenter);
      break;
    case 7:
      /* case 7, Q = sum_i S  U  cos(q-q*) */
      geneF = -1.0 * bias->cvstrength * deltaCV * bias->signU * bias->bic[indexbic].cvcoeff
              * sin(bias->bic[indexbic].val - bias->bic[indexbic].cvcenter);
      break;
    default:
      fprintf(stderr, "\nHL Error: calc_gene_biasForce unknown gbf %d\n\n", gbf);
      exit(1);
  }

  /* the do_ddih_HL() calculate -1.0 * partial_q over partial_x */
  /* thus we need to flip the sign of generalized bias force as well */
  /* the force calculated above are negative of generalized bias force */

  return geneF;
}

/* B-matrix elements: partial_dihedral over partial_cartesianCoordinats */
/* do_ddih_HL() give negative sign B-matrix elements */
void do_ddih_HL(const rvec xi,const rvec xj,const rvec xk,const rvec xl, rvec ddih[4])
{
  /* this is NEGATIVE partial_dih over partial_x */
  rvec r_ij, r_kj, r_kl;
  rvec uvec,vvec,svec, m,n;
  real iprm,iprn,nrkj,nrkj2;
  real a,p,q;
  int i, j;

  for(i=0; i<4; i++) {
    for(j=0; j<3; j++) {
      ddih[i][j] = 0;
    }
  }
  rvec_sub(xi, xj, r_ij);
  rvec_sub(xk, xj, r_kj);
  rvec_sub(xk, xl, r_kl);

  cprod(r_ij,r_kj,m);                   /*  9           */
  cprod(r_kj,r_kl,n);

  iprm  = iprod(m,m);           /*  5   */
  iprn  = iprod(n,n);           /*  5   */
  nrkj2 = iprod(r_kj,r_kj);     /*  5   */
//  nrkj  = nrkj2*gmx_invsqrt(nrkj2);     /* 10   */
  nrkj  = nrkj2*(gmx::invsqrt(nrkj2));     /* 10   */

  a     = -nrkj/iprm;   /* 11   */
  svmul(a,m,ddih[0]);           /*  3   */
  a     = nrkj/iprn;    /* 11   */
  svmul(a,n,ddih[3]);           /*  3   */
  p     = iprod(r_ij,r_kj);     /*  5   */
  p    /= nrkj2;                /* 10   */
  q     = iprod(r_kl,r_kj);     /*  5   */
  q    /= nrkj2;                /* 10   */
  svmul(p,ddih[0],uvec);                /*  3   */
  svmul(q,ddih[3],vvec);                /*  3   */
  rvec_sub(uvec,vvec,svec);     /*  3   */
  rvec_sub(svec, ddih[0], ddih[1]);     /*  3   */
  rvec_add(svec, ddih[3], ddih[2]);     /*  3   */
  svmul(-1, ddih[2], ddih[2]);  /*  3   */
}

/* For one dihedral, calculate bias force on Cartesian coordinate, given generalized force */
void calc_C_biasForce(BIC *bic, rvec *biasf, rvec *x)
{
  double geneF;
  int ai, aj, ak, al;
  rvec ddih[4];
  int i, m, atomid;

  geneF = bic->neggenef;
  ai = bic->atom[0];
  aj = bic->atom[1];
  ak = bic->atom[2];
  al = bic->atom[3];

  /* calculate B-matrix elements */
  /* do_ddih_HL() give negative sign for B-matrix elements */
  do_ddih_HL(x[ai], x[aj], x[ak], x[al], ddih);

  for (i = 0; i < 4; i++){
    atomid = bic->atom[i];
    for (m = 0; m < 3; m++){
      biasf[atomid][m] += (geneF * ddih[i][m]);
    }
  }
}

/* add bias force on collective variable */
double add_bias_CV(BIAS *bias, rvec *f, rvec *x, double t, int hlstep)
{
  double biasPot;
  double dq;
  int i;
 
  /* calculate value of collective variable */
  bias->cv = calcCV(bias->numic, bias->bic, bias->gbf, bias->signU); 
  /* use first step collective variable value and target value */
  /* to determine bias force direction on CV */
  if (hlstep == 0){
    if (bias->cv < bias->cvtarget){
      bias->cvdirection = 1;
      bias->cvfastbin = 0;
    } else{
      bias->cvdirection = -1;
      bias->cvfastbin = bias->cvnumbin - 1;
    }
  }
  
  /* Do not apply bias force, if system ever reach production basin */
  if (bias->everbasin == 1){
    fprintf(stderr, "\nHL no bias force on CV due to everbasin == 1, at t = %lf ps\n\n", t);
    return 0.0;       
  }
 
  /* CV bias bin, which CV located at */
  bias->cvcurrbin = whichbiasbin(bias->cv, bias->cvnumbin, bias->cvbd);
  if ((bias->cvcurrbin == -100) || (bias->cvcurrbin == -1)){
    fprintf(stderr, "\nHL error: add_bias_CV, %lf out of CV range [%lf %lf]\n\n",
            bias->cv, bias->cvbd[0], bias->cvbd[bias->cvnumbin]);
    exit(1);
  }
 
  if (bias->cvdirection == 1){
    bias->cveq = bias->cvbd[bias->cvcurrbin+1];
  } else if (bias->cvdirection == -1){
    bias->cveq = bias->cvbd[bias->cvcurrbin];
  } else{
    fprintf(stderr, "\nHL error: add_bias_CV, wrong cvdirection %5d\n\n", bias->cvdirection);
    exit(1);
  }

  /* we use lrbb for control using left or right boundary of bias bin to */
  if ((bias->gbf == 5) || (bias->gbf == 6) || (bias->gbf == 7)){
    if (bias->lrbb == 1){
      /* use right boundary */
      bias->cveq = bias->cvbd[bias->cvcurrbin+1];
    } else if (bias->lrbb == -1){
      /* use left boundary */
      bias->cveq = bias->cvbd[bias->cvcurrbin];
    } else{
      fprintf(stderr, "\nHL error: gbf lrbb %5d %5d\n\n", bias->gbf, bias->lrbb);
      exit(1);
    }
  } 

  /* calculate bias potential and generalized bias force */
  dq = bias->cv - bias->cveq;
  biasPot = 0.5 * bias->cvstrength * dq * dq + 
            (bias->cvnumbin - 1 - bias->cvcurrbin) * 0.5 * bias->cvstrength *
            bias->cvbinsize * bias->cvbinsize;

  /* for each dihedral, calculate generalized force on it */
  /* then use negative B-matrix elements from do_ddih_HL() to get bias force on Cartesian coordinate */
  for (i = 0; i < bias->numic; i++){
    bias->bic[i].neggenef = calc_gene_biasForce(bias, dq, bias->gbf, i); 
    calc_C_biasForce(&(bias->bic[i]), bias->biasf, x);
  }

  return biasPot;
}

/* bias force are saved in bias->bias, which is vector of double  */
/* add this bias force to rvec *f */
void bias_to_rvecf(rvec *biasf, int natom, rvec *f)
{
  int i, m;

  for (i = 0; i < natom; i++){
    for (m = 0; m < 3; m++){
      f[i][m] += biasf[i][m];
    }
  } 
}

/* output bias force */
void outputBiasForce(BIAS *bias, double t)
{
  int i;
  fwrite(&t, sizeof(double), 1, bias->fpbf);
  for (i = 0; i < bias->natom; i++){
    fwrite(bias->biasf[i], sizeof(double), 3, bias->fpbf);
  }
}

/* output collective variable information */
void outputCV(BIAS *bias, double t)
{
  //fprintf(bias->fpbias, "%lf %16.8f %16.8f %5d ", t, bias->cv, bias->cveq, bias->cvcurrbin);
  //fprintf(bias->fpbias, "\n");
}

/* output whether system reachs target basin, value of order parameters */
void outputBasin(BIAS *bias, double t)
{
  int i;

  fprintf(bias->fpbasin, "%lf ", t);
  for (i = 0; i < bias->numic; i++){
    fprintf(bias->fpbasin, "%12.8e ", bias->bic[i].val);
  }
  for (i = 0; i < bias->numic; i++){
    fprintf(bias->fpbasin, "%12.8e ", -1.0 * bias->bic[i].neggenef);
  }
  fprintf(bias->fpbasin, "\n");

/*
  fprintf(bias->fpbasin, "%lf %5d %5d ", t, bias->whereach, bias->everbasin);
  for (i = 0; i < bias->numic; i++){
    if (bias->bic[i].od == 1){
      fprintf(bias->fpbasin, "%lf ", bias->bic[i].val);
    }
  }
  fprintf(bias->fpbasin, "\n");
 */
}

/* output bias information */
void outputBIAS(BIAS *bias, double t, int hlstep)
{
  if ((hlstep % (bias->deffreq)) == 0){
    //outputBiasForce(bias, t);
    //outputBasin(bias, t);
  } 

  switch (bias->whecv){
    case 1:
      outputCV(bias, t);
      break;
    default:
      fprintf(stderr, "\nHL Error: outputBIAS unknown whecv %d\n\n", bias->whecv);
      exit(1);
  }
}

/* whether system reaches target basin, according to order parameters */
/* 1 for yes, 0 for no */
int reachbasin(BIAS *bias)
{
  int i;

  /* as long as one order parameter is not in target range, then system is  not in target basin */
  for (i = 0; i < bias->numic; i++){
    if (bias->bic[i].od == 1){  /* this dihedral is order parameter */
      if ((bias->bic[i].val < bias->bic[i].odlow) || (bias->bic[i].val > bias->bic[i].odup)){
        return 0;
      }
    }
  }

  /* if all order parameters locate within target range, then system reaches target basin */
  return 1;
}

void calcU(int numdih, int **atom, double *dih, rvec *x, FILE *fp, double t, rvec boxsize)
{
  int ai, aj, ak, al;
  int i;

  fwrite(&t, sizeof(double), 1, fp);
  for (i = 0; i < numdih; i++){
    ai = atom[i][0];
    aj = atom[i][1];
    ak = atom[i][2];
    al = atom[i][3];
    dih[i] = dih_angle_HL(x[ai], x[aj], x[ak], x[al], boxsize);
  }
  fwrite(dih, sizeof(double), numdih, fp);
}

/* the root function for adding bias force */
void add_bias_force(t_inputrec *ir, rvec *f, rvec *x, double t, int hlstep)
{
  BIAS *bias;
  double biasPot;
  int i; 

  bias = &(ir->bias);

  /* calculate value of BBdih involved in U coordiantes */
  if (bias->wheU > 0){
    calcU(bias->wheU, bias->uatom, bias->ubbdih, x, bias->fpdih, t, bias->boxsize);
  }

  /* calculate value of each dihedral first, this is required by all different bias methods */
  for (i = 0; i < bias->numic; i++){
    /* fprintf(stderr, "IC %d ", i); * testing purpose only **/
    calcdihBIC(x, &(bias->bic[i]), bias->boxsize);
  }

  /* set bias->biasf to zero */
  setzeroHL(bias->biasf, bias->natom);
  
  /* calculate bias force with different analytical expression, according to bias->bgf */
  /* bias force on Cqrtesian coordinates are saved in bias->biasf */
  switch(bias->whecv) {
    case 1:
      /* bias force on collective variable */    
      biasPot = add_bias_CV(bias, f, x, t, hlstep);
      break;
    default:
      fprintf(stderr, "\nError: no implimentation for whecv %d\n\n", bias->whecv);
  }

  /* add bias force (bias->biasf) to the rvec *f */
  bias_to_rvecf(bias->biasf, bias->natom, f);

  /* according to order parameters, check whether system reach bias basin */
  bias->whereach = reachbasin(bias); 
  if (bias->whereach == 1){
    bias->everbasin = 1;
  }

  /* whether stop simulation when system reaches target basin */
  if ((bias->whebiasstop == 1) && (bias->whereach == 1)){
    bias->bLastStep = TRUE;
  }

  /* output */
  outputBIAS(bias, t, hlstep);
}

/* end of adding by Huiyu Li */

