### environment

import scanpy as sc
import scvi
# import scarches as sca

### global settings

# UserWarning: Since v1.0.0, scvi-tools no longer uses a random seed by default
scvi.settings.seed = 0
TESTING_MODE = True

### functions

def run_scVI_reference(adata, out_mod): 
    """
    Learns a latent representation space and saves it into adata.obsm
    
    Parameters
    ----------
    adata
        input adata object 
    out_mod
        output file for scVI model   
    """
    # NotImplementedError: scArches currently does not support models with extra categorical covariates.
    scvi.model.SCVI.setup_anndata(
        adata,
        layer = "counts"
    )

    vae = scvi.model.SCVI(
        adata, 
        use_layer_norm = "both",
        use_batch_norm = "none",
        encode_covariates = True,
        dropout_rate = 0.2,
        n_layers = 2)
    if TESTING_MODE:
        vae.train(max_epochs = 2) # FOR TESTING
    else:
        vae.train()
    vae.save(out_mod, overwrite = True)


def run_scANVI_reference(adata, in_mod, out_mod, meta):
    """
    Learns a latent representation space and saves it into adata.obsm
    
    Parameters
    ----------
    adata
        input adata object
    in_mod
        input file for scVI model      
    out_mod
        output file for scANVI model 
    meta
        name of the observation covariate  
    """    
    vae_ref = scvi.model.SCVI.load(
        in_mod, 
        adata
    )

    adata.obs["labels_scanvi"] = adata.obs[meta].values

    vae_ref_scan = scvi.model.SCANVI.from_scvi_model(
        vae_ref, 
        unlabeled_category = "Unknown",
        labels_key = "labels_scanvi"
    )
    max_epochs = 20
    if TESTING_MODE:
        max_epochs = 2 # FOR TESTING
    vae_ref_scan.train(
        max_epochs = max_epochs,
        n_samples_per_label = 100
    )
    vae_ref_scan.save(out_mod, overwrite = True)


def run_scVI_query(adata, ref_mod, q_mod, ref_name):
    """
    Learns a latent representation space and saves it
    
    Parameters
    ----------
    adata
        input adata query object
    ref_mod
        input file for reference scVI model      
    q_mod
        output file for query scVI model 
    ref_name
        reference identifier
    """   
    vae_q = scvi.model.SCVI.load_query_data(
        adata,
        ref_mod
    )
    max_epochs = 20
    if TESTING_MODE:
        max_epochs = 2 # FOR TESTING
    vae_q.train(
        max_epochs = max_epochs,
        plan_kwargs = {"weight_decay": 0.0}, 
        check_val_every_n_epoch = 10
    ) 
    vae_q.save(q_mod, overwrite = True)


def run_scANVI_query(adata, ref_mod, q_mod, meta, ref_name):
    """
    Learns a latent representation space and saves it
    
    Parameters
    ----------
    adata
        input adata query object
    ref_mod
        input file for reference scANVI model      
    q_mod
        output file for query scANVI model 
    meta
        name of the observation covariate
    ref_name
        reference identifier
    """       
    adata.obs["labels_scanvi"] = "Unknown"

    scvi.model.SCANVI.prepare_query_anndata(
        adata, 
        ref_mod
    ) 
    vae_q = scvi.model.SCANVI.load_query_data(
        adata,
        ref_mod
    )
    max_epochs = 100
    if TESTING_MODE:
        max_epochs = 10 # FOR TESTING
    vae_q.train(
        max_epochs = 100,
        plan_kwargs = {"weight_decay": 0.0}, 
        check_val_every_n_epoch = 10
    ) 
    vae_q.save(q_mod, overwrite = True)


# separate the model calculation and the adata update for all the other cases
# the rationale is that the model will take a lot of time to compute and writing of adata to disk cannot be done concurrently

def update_adata_scVI_reference(adata, ref_mod):
    """
    Saves the model in adata
    
    Parameters
    ----------
    adata
        input reference adata object
    ref_mod
        input reference scVI model      
    """   
    model = scvi.model.SCVI.load(ref_mod, adata)
    adata.obsm["X_scvi"] = model.get_latent_representation()
    adata.layers["normalized_scvi"] = model.get_normalized_expression(library_size=10e6)
    return adata


def update_adata_scANVI_reference(adata, ref_mod, meta):
    """
    Saves the model in adata
    
    Parameters
    ----------
    adata
        input reference adata object
    ref_mod
        input reference scANVI model
    meta
        name of the observation covariate
    """
    model = scvi.model.SCANVI.load(ref_mod, adata)
    latent_id = "X_scanvi_" + meta
    adata.obsm[latent_id] = model.get_latent_representation()
    return adata


def update_adata_scVI_query(adata, q_mod, ref_name):
    """
    Saves the model in adata
    
    Parameters
    ----------
    adata
        input query adata object     
    q_mod
        input query scANVI model 
    ref_name
        reference identifier
    """   
    model = scvi.model.SCVI.load(q_mod, adata)
    latent_id = "X_scvi_" + ref_name
    adata.obsm[latent_id] = model.get_latent_representation()
    return adata


def update_adata_scANVI_query(adata, q_mod, ref_name, meta):
    """
    Saves the model in adata
    
    Parameters
    ----------
    adata
        input query adata object  
    q_mod
        input query scANVI model 
    ref_name
        reference identifier
    meta
        name of the observation covariate
    """   
    model = scvi.model.SCANVI.load(q_mod, adata)
    latent_id = "X_scanvi_" + ref_name + "_" + meta
    adata.obsm[latent_id] = model.get_latent_representation()
    return adata


def predict_scANVI_query_labels(adata, q_mod):
    """
    Predict cell type labels
    
    Parameters
    ----------
    adata
        input query adata object  
    q_mod
        input query scANVI model 
    """ 
    model = scvi.model.SCANVI.load(q_mod, adata)
    labels_hard = model.predict()
    labels_soft = model.predict(soft = True)
    return labels_hard, labels_soft


