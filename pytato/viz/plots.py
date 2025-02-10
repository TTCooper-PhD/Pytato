import matplotlib.pyplot as plt
import seaborn as sns
import logomaker
import pandas as pd

class Plots:
    """
    Visualization methods for peptide/protein data
    """
    @staticmethod
    def sequence_logo(data, width=0.9, vpad=0.2, stack_order='small_on_top', 
                     color_scheme='hydrophobicity', font_name="Ebrima", 
                     figsize=[6,2], **kwargs):
        """
        Create a sequence logo plot
        """
        logo = logomaker.Logo(data,
                            width=width,
                            vpad=vpad,
                            stack_order=stack_order,
                            color_scheme=color_scheme,
                            font_name=font_name,
                            fade_probabilities=False,
                            center_values=False,
                            alpha=1,
                            figsize=figsize)
                            
        # Style using Logo methods
        logo.style_spines(spines=['left', 'right'], visible=True)
        
        # Style using Axes methods
        logo.ax.set_xticks(range(len(data)))
        logo.ax.set_xticklabels('%+d'%x for x in range(-5, 6))
        logo.ax.set_yticks([0, .5, 1])
        logo.ax.axvline(4.5, color='k', linewidth=2, linestyle=':')
        logo.ax.set_ylabel('probability', fontsize=12)
        
        return logo 
    
    @staticmethod

    def peptide_hex(data,char1="m/z",char2="GRAVY",figsize=(5,4),color="#4CB391"):
        plt.figure(figsize=figsize)
        x=data[char1]
        y=data[char2]
        plt.figure(figsize=figsize)
        sns.jointplot(x=x,y=y,kind="hex",color=color)
        plt.show()

    @staticmethod
    def double_kde(dataframe1,dataframe2, selected_column='z', label1='Unmodified', label2='Modified',common_norm=True,figsize=(8,5)):
        plt.figure(figsize=figsize)
        sns.kdeplot(data=dataframe1, x=selected_column, label=label1, common_norm=common_norm)
        sns.kdeplot(data=dataframe2, x=selected_column, label=label2, common_norm=common_norm)
        plt.xlabel(selected_column)
        plt.ylabel('Density')
        plt.title(f'Distribution of {selected_column} Values: {label1} vs {label2}')
        plt.legend()
        plt.show()
