# GOTaxonConstraint
The Taxon Constraints for Gene Onotology Terms

We have implemented a method of representing taxon constraints in terms of evolutionary gain and loss events relative to the species tree of life, and inferred constraints for the majority of GO by combining manual curation and rule based propagation.


We first manually curate a seed list of high levels GO terms based on biological knowledge, making use of a broad range of ontologies that are linked to the GO and have not previously been applied to this problem. We propagate these constraints to related more specific GO terms based on the “true path rule” of the GO hierarchy. The constraints are then compared with experimental annotations to find conflicts. We manually check the conflicts and the GO terms without taxon constraints. Then we modify the seed list and begin another round of propagation. This process was repeated for several cycles of curation, propagation and correction. 


By expressing taxon constraints in terms of evolutionary gains and losses, we can apply the existing knowledge about evolutionary histories, such as the evolution of eukaryotic cell, multicellularity and lineage-specific elaboration of anatomical structures, to the problem of constructing a computational representation of biological systems. We hope these constraints will provide valuable information for both GO curators and the wider community of GO and GO annotation users, and will prove useful in improving the accuracy of both manual and, in particular, computationally predicted, GO annotations.

If you use the GO taxon constraints, please cite our paper
GOTaxon: Representing the evolution of biological functions in the Gene Ontology.


******************************************************
The GO taxon constraints are included in file GOTaxonConstraints.txt
******************************************************
The format is like this:

           GO id => ">Gain|" followd by one or more NCBITaxon ids separated by ";" which mean this GO term is constrained to one of these GO terms (relationship is "Or")
           Some GO terms also have "Loss", which is represented by ">Loss|" followed by one or more NCBITaxon ids.

          example: 'GO:0000578' => '>Gain|NCBITaxon:3193;NCBITaxon:33213;NCBITaxon:5782;NCBITaxon:451864;>Loss|NCBITaxon:4896;NCBITaxon:4892;',


          There are also some GO terms that have representation: ">Chebi|" followed by NCBITaxon ids. These representation indicates that the taxon constraint for this term comes from Chebi info only.

           example: 'GO:0000718' => '>Chebi|NCBITaxon:1',


******************************************************
For simplicity of usage, the GO taxon constraints are also presented in a differnt format.
In file GOTaxonConstraintsExtantSpecies.csv
The file indicates whether a GO term should be annotated to certain extant species. 
******************************************************

The format is like this:
         each line lists if the GO term in the first column could be annotated for each species in other columns 
         The first row lists extant species in 5 letter syns, for details check http://pantherdb.org/panther/summaryStats.jsp and https://www.uniprot.org/help/taxonomy)				
         
         example:
         ## For the numbers in the table:																	
         ## 2 means a GO term could be annotated to the species
         ## -1 means a GO term should not be annotated to the species				
         ## 0 also means a GO term should not be annotated to the species (Different codes used as reasons for forbidden annotation are different. You could treat 0 and -1 as the same)
         
         GOid	       CANAL	TAKRU	ORYSJ
         GO:1903097		2	2	2
         GO:2001001		2	2	2
         GO:0006285		2	2	2
         GO:0015777		2	2	2
         GO:0061408		-1	-1	-1

This means GO:2001001 could be annotated with CANAL(Candida albicans), TAKRU (Takifugu rubripes) and ORYCJ (Oryza sativa subsp. japonica (Rice)),  while GO:0061108 could not be annotated with these species.


