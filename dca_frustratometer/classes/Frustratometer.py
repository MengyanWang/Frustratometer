"""Provide the primary functions."""
from pathlib import Path
from .. import pdb
from .. import dca

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

#Import other modules
from ..utils import _path
from .. import frustration

__all__=['PottsModel']
##################
# PFAM functions #
##################


# Class wrapper
class Frustratometer:

    # @property
    # def sequence_cutoff(self):
    #     return self._sequence_cutoff

    # @sequence_cutoff.setter
    # def sequence_cutoff(self, value):
    #     self.mask = frustration.compute_mask(self.distance_matrix, self.distance_cutoff, self.sequence_cutoff)
    #     self._sequence_cutoff = value
    #     self._native_energy = None
    #     self._decoy_fluctuation = {}

    # @property
    # def distance_cutoff(self):
    #     return self._distance_cutoff

    # @distance_cutoff.setter
    # def distance_cutoff(self, value):
    #     self.mask = frustration.compute_mask(self.distance_matrix, self.distance_cutoff, self.sequence_cutoff)
    #     self._distance_cutoff = value
    #     self._native_energy = None
    #     self._decoy_fluctuation = {}

    def native_energy(self,sequence:str = None,ignore_contacts_with_gaps:bool = False):
        if sequence is None:
            sequence=self.sequence
        else:
            self._native_energy=frustration.compute_native_energy(sequence, self.potts_model, self.mask,ignore_contacts_with_gaps)
        if not self._native_energy:
            self._native_energy=frustration.compute_native_energy(sequence, self.potts_model, self.mask,ignore_contacts_with_gaps)
        return self._native_energy

    def sequences_energies(self, sequences:np.array, split_couplings_and_fields:bool = False):
        return frustration.compute_sequences_energy(sequences, self.potts_model, self.mask, split_couplings_and_fields)

    def fields_energy(self, sequence=None):
        if sequence is None:
            sequence=self.sequence
        return frustration.compute_fields_energy(sequence, self.potts_model)

    def couplings_energy(self, sequence:str = None,ignore_contacts_with_gaps:bool = False):
        if sequence is None:
            sequence=self.sequence
        return frustration.compute_couplings_energy(sequence, self.potts_model, self.mask,ignore_contacts_with_gaps)
        
    def decoy_fluctuation(self, sequence:str = None,kind:str = 'singleresidue',mask:np.array = None):
        if sequence is None:
            sequence=self.sequence
            if kind in self._decoy_fluctuation:
                return self._decoy_fluctuation[kind]
        if not isinstance(mask, np.ndarray):
            mask=self.mask
        if kind == 'singleresidue':
            fluctuation = frustration.compute_singleresidue_decoy_energy_fluctuation(sequence, self.potts_model, mask)
        elif kind == 'mutational':
            fluctuation = frustration.compute_mutational_decoy_energy_fluctuation(sequence, self.potts_model, mask)
        elif kind == 'configurational':
            fluctuation = frustration.compute_configurational_decoy_energy_fluctuation(sequence, self.potts_model, mask)
        elif kind == 'contact':
            fluctuation = frustration.compute_contact_decoy_energy_fluctuation(sequence, self.potts_model, mask)
        else:
            raise Exception("Wrong kind of decoy generation selected")
        self._decoy_fluctuation[kind] = fluctuation
        return self._decoy_fluctuation[kind]

    def decoy_energy(self, kind:str = 'singleresidue'):
        return self.native_energy() + self.decoy_fluctuation(kind=kind)

    def scores(self):
        return frustration.compute_scores(self.potts_model)

    def frustration(self, sequence:str = None, kind:str = 'singleresidue', mask:np.array = None, aa_freq:np.array = None, correction:int = 0):
        if sequence is None:
            sequence=self.sequence
        if not isinstance(mask, np.ndarray):
            mask=self.mask
        decoy_fluctuation = self.decoy_fluctuation(sequence=sequence,kind=kind, mask=mask)
        if kind == 'singleresidue':
            if aa_freq is not None:
                aa_freq = self.aa_freq
            return frustration.compute_single_frustration(decoy_fluctuation, aa_freq, correction)
        elif kind in ['mutational', 'configurational', 'contact']:
            if aa_freq is not None:
                aa_freq = self.contact_freq
            return frustration.compute_pair_frustration(decoy_fluctuation, aa_freq, correction)

    def plot_decoy_energy(self, kind:str = 'singleresidue', method:str = 'clustermap'):
        native_energy = self.native_energy()
        decoy_energy = self.decoy_energy(kind=kind)
        if kind == 'singleresidue':
            g = frustration.plot_singleresidue_decoy_energy(decoy_energy, native_energy, method)
            return g

    def roc(self):
        return frustration.compute_roc(self.scores(), self.distance_matrix, self.distance_cutoff)

    def plot_roc(self):
        frustration.plot_roc(self.roc())

    def auc(self):
        """Computes area under the curve of the receiver-operating characteristic.
           Function intended"""
        return frustration.compute_auc(self.roc())

    def vmd(self, single:str = 'singleresidue', pair:str = 'mutational', aa_freq:np.array = None, correction:int = 0, max_connections:int = 100):
        tcl_script = frustration.write_tcl_script(self.pdb_file, self.chain,
                                      self.frustration(single, aa_freq=aa_freq, correction=correction),
                                      self.frustration(pair, aa_freq=aa_freq, correction=correction),
                                      max_connections=max_connections)
        frustration.call_vmd(self.pdb_file, tcl_script)

    def view_frustration(self, single:str = 'singleresidue', pair:str = 'mutational', aa_freq:np.array = None, correction:int = 0, max_connections:int = 100):
        import py3Dmol
        pdb_filename = self.pdb_file
        shift=self.init_index_shift+1
        pair_frustration=self.frustration(kind=pair)*np.triu(self.mask)
        residues=np.arange(len(self.sequence))
        r1, r2 = np.meshgrid(residues, residues, indexing='ij')
        sel_frustration = np.array([r1.ravel(), r2.ravel(), pair_frustration.ravel()]).T
        minimally_frustrated = sel_frustration[sel_frustration[:, -1] > 1]
        frustrated = sel_frustration[sel_frustration[:, -1] < -.78]
        
        view = py3Dmol.view(js='https://3dmol.org/build/3Dmol.js')
        view.addModel(open(pdb_filename,'r').read(),'pdb')

        view.setBackgroundColor('white')
        view.setStyle({'cartoon':{'color':'white'}})
        
        for i,j,f in frustrated:
            view.addLine({'start':{'chain':'A','resi':[str(i+shift)]},'end':{'chain':'A','resi':[str(j+shift)]},
                        'color':'red', 'dashed':False,'linewidth':3})
        
        for i,j,f in minimally_frustrated:
            view.addLine({'start':{'chain':'A','resi':[str(i+shift)]},'end':{'chain':'A','resi':[str(j+shift)]},
                        'color':'green', 'dashed':False,'linewidth':3})

        view.zoomTo(viewer=(0,0))

        return view

    def view_frustration_pair_distribution(self,kind:str ="singleresidue",include_long_range_contacts:bool =True):
        #Ferrerio et al. (2007) pair distribution analysis included long-range contacts.
        if include_long_range_contacts==True:
            mask=frustration.compute_mask(self.distance_matrix, distance_cutoff=None, sequence_distance_cutoff = 2)
        else:
            mask=self.mask

        frustration_values=self.frustration(kind=kind,mask=mask)

        residue_ca_coordinates=(self.structure.select('calpha').getCoords())
        all_pairs=[]
        for i in range(0,frustration_values.shape[0]):
            for j in range(0,frustration_values.shape[0]): 
                all_pairs.append([i,j]) 
        all_pairs=np.array(all_pairs)
        print(all_pairs)
        if kind=="singleresidue":
            frustration_values = np.array(frustration_values)*4
        elif kind in ['mutational', 'configurational']:
            frustration_values = frustration_values.flatten()*4
        print(frustration_values)
        minimally_fustrated = frustration_values>0.78
        neutrally_fustrated = (frustration_values>=-1) & (frustration_values<=0.78) 
        maximally_fustrated = frustration_values<-1

        pair_groups={'All':all_pairs,
                    'Minimally frustrated':all_pairs[minimally_fustrated],
                    'Neutrally frustrated':all_pairs[neutrally_fustrated],
                    'Maximally frustrated':all_pairs[maximally_fustrated],}  
        print(pair_groups)
        ###
        for group,pairs in pair_groups.items():
            # Calculate centers
            centers = (residue_ca_coordinates[pairs[:, 0]] + residue_ca_coordinates[pairs[:, 1]]) # Px3 array
            if kind in ['mutational', 'configurational']:
                centers=centers/2
            #Calculate distances between centers
            diff = centers[:, np.newaxis] - centers 
            distances = np.sqrt((diff ** 2).sum(axis=-1)) # PxP array

            # Set a step size (dr) and a maximum radius (max_r)
            n_bins = 100
            max_r = 30.0

            #Number of points
            N = len(pairs) # if needed

            # Calculate the volume of a spherical shell of thickness dr
            r=np.linspace(0,max_r,n_bins)
            r_m=(r[1:]+r[:-1])/2
            shell_vol = 4/3 * np.pi * (r[1:]**3-r[:-1]**3)

            # Histogram the distances; these are our "raw" g(r)
            dist_hist = np.histogram(distances, bins=r)[0]
            # Normalize the histogram
            rdf_pairs = dist_hist / (N)
            if kind in ['mutational', 'configurational']:
                rdf_pairs=rdf_pairs/shell_vol
            plt.plot(r_m,rdf_pairs,label=group)

        plt.xlim(0.01,20)
        plt.legend() 


    def view_frustration_histogram(self,kind:str = "singleresidue"):
        def frustration_type(x):
            if x>.78:
                frustration_class="Minimally Frustrated"
            elif x<-1:
                frustration_class="Frustrated"
            else:
                frustration_class="Neutral"
            return frustration_class

        frustration_values=self.frustration(kind=kind)
        cb_distance_matrix=self.distance_matrix
        i_index,j_index=np.triu_indices(cb_distance_matrix.shape[0], k = 1)

        flattened_cb_distance_matrix=cb_distance_matrix[i_index,j_index]
        contact_pairs=list((zip(i_index,j_index)))
        ###
        if kind=="singleresidue":
            frustration_dataframe=pd.DataFrame(data=np.array([range(0,len(frustration_values)),frustration_values]).T,
                                    columns=["Residue Index","F_i"])
            #Classify frustration type
            frustration_dataframe["F_i_Type"]=frustration_dataframe["F_i"].apply(lambda x: frustration_type(x))

            #Plot histogram of all frustration values.
            with sns.plotting_context("poster"):
                plt.figure(figsize=(10,5))

                g=sns.histplot(data=frustration_dataframe,x="F_i",bins=20)
                ymin, ymax = g.get_ylim()
                g.vlines(x=[-1, .78], ymin=ymin, ymax=ymax, colors=['gray', 'gray'], ls='--', lw=2)
                plt.title(f"N={len(frustration_dataframe)}")
                plt.show()
            ###
            print(f"{((len(frustration_dataframe.loc[frustration_dataframe['F_i_Type']=='Minimally Frustrated'])/len(frustration_dataframe))*100):.2f}% Residues are Minimally Frustrated")
            print(f"{((len(frustration_dataframe.loc[frustration_dataframe['F_i_Type']=='Frustrated'])/len(frustration_dataframe))*100):.2f}% Residues are Frustrated")
            print(f"{((len(frustration_dataframe.loc[frustration_dataframe['F_i_Type']=='Neutral'])/len(frustration_dataframe))*100):.2f}% Residues are Neutral")

        elif kind in ['mutational', 'configurational']:
            frustration_values=frustration_values[i_index,j_index]
            
            frustration_dataframe=pd.DataFrame(data=np.array([contact_pairs,flattened_cb_distance_matrix,frustration_values]).T,
                                                columns=["i,j Pair","Original_Distance_ij","F_ij"])
            frustration_dataframe=frustration_dataframe.dropna(subset=["F_ij"])
            #Classify frustration type
            frustration_dataframe["F_ij_Type"]=frustration_dataframe["F_ij"].apply(lambda x: frustration_type(x))
            frustration_dataframe["Contact_Type"]=np.where(frustration_dataframe["Original_Distance_ij"]<6.5,"Direct","Water-Mediated")
            #Plot histogram of all frustration values.
            with sns.plotting_context("poster"):
                fig,axes=plt.subplots(1,2,figsize=(15,5),sharex=True)

                g=sns.histplot(data=frustration_dataframe.loc[frustration_dataframe["Contact_Type"]=="Direct"],x="F_ij",bins=20,ax=axes[0])
                ymin, ymax = g.get_ylim()
                g.vlines(x=[-1, .78], ymin=ymin, ymax=ymax, colors=['gray', 'gray'], ls='--', lw=2)
                axes[0].title.set_text(f"Direct Contacts (N={len(frustration_dataframe.loc[frustration_dataframe['Contact_Type']=='Direct'])})")
                ###
                g=sns.histplot(data=frustration_dataframe.loc[frustration_dataframe["Contact_Type"]=="Water-Mediated"],x="F_ij",bins=20,ax=axes[1])
                ymin, ymax = g.get_ylim()
                g.vlines(x=[-1, .78], ymin=ymin, ymax=ymax, colors=['gray', 'gray'], ls='--', lw=2)
                axes[1].title.set_text(f"Protein- & Water-\nMediated Contacts (N={len(frustration_dataframe.loc[frustration_dataframe['Contact_Type']=='Water-Mediated'])})")

                plt.tight_layout()
                plt.show()
            ###
            direct_contact_frustration_dataframe=frustration_dataframe.loc[frustration_dataframe["Contact_Type"]=="Direct"]
            print(f"{((len(direct_contact_frustration_dataframe.loc[direct_contact_frustration_dataframe['F_ij_Type']=='Minimally Frustrated'])/len(direct_contact_frustration_dataframe))*100):.2f}% Direct Contacts are Minimally Frustrated")
            print(f"{((len(direct_contact_frustration_dataframe.loc[direct_contact_frustration_dataframe['F_ij_Type']=='Frustrated'])/len(direct_contact_frustration_dataframe))*100):.2f}% Direct Contacts are Frustrated")
            print(f"{((len(direct_contact_frustration_dataframe.loc[direct_contact_frustration_dataframe['F_ij_Type']=='Neutral'])/len(direct_contact_frustration_dataframe))*100):.2f}% Direct Contacts are Neutral")





