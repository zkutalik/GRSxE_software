get_fY_results  =  function( rslurm_jobname = NULL,
                             fY_results     = NULL ){
    if (is.null( fY_results )) {
        require( rslurm )
        fY_results  =  get_slurm_out( rslurm_jobname )
    }

}














