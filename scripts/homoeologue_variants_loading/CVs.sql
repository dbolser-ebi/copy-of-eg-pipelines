-- atrib_type table

INSERT INTO `attrib_type` VALUES (450,'A_homoeologue_locus','homoeologous locus in A','homoeologous locus in A');
INSERT INTO `attrib_type` VALUES (451,'B_homoeologue_locus','homoeologous locus in B','homoeologous locus in B');
INSERT INTO `attrib_type` VALUES (452,'D_homoeologue_locus','homoeologous locus in D','homoeologous locus in D');

-- attrib table

INSERT INTO `attrib` VALUES (1004,9,'Inter-homoeologous variants between A and B, not present in D');
INSERT INTO `attrib` VALUES (1005,9,'Inter-homoeologous variants between A and D, not present in B');
INSERT INTO `attrib` VALUES (1006,9,'Inter-homoeologous variants between B and D, not present in A');
INSERT INTO `attrib` VALUES (1007,9,"Inter-homoeologous variants between A, B and D, where A, B and D don't share the same allele");
INSERT INTO `attrib` VALUES (1008,9,"Inter-homoeologous variants between A, B and D, where A and B share the same allele");
INSERT INTO `attrib` VALUES (1009,9,"Inter-homoeologous variants between A, B and D, where A and D share the same allele");
INSERT INTO `attrib` VALUES (1010,9,"Inter-homoeologous variants between A, B and D, where B and D share the same allele");
INSERT INTO `attrib` VALUES (1011,9,"All SNPs from CerealsDB genotyping platforms");
INSERT INTO `attrib` VALUES (1012,9,"All variants classified as inter-homoeologues");

-- variation_set table

INSERT INTO `variation_set` VALUES (4,'Inter-homoeologous variants between A and B, not present in D','Inter-homoeologous variants between A and B bread wheat component genomes, not present in D',1004);
INSERT INTO `variation_set` VALUES (5,'Inter-homoeologous variants between A and D, not present in B','Inter-homoeologous variants between A and D bread wheat component genomes, not present in B',1005);
INSERT INTO `variation_set` VALUES (6,'Inter-homoeologous variants between B and D, not present in A','Inter-homoeologous variants between B and D bread wheat component genomes, not present in A',1006);
INSERT INTO `variation_set` VALUES (7,"Inter-homoeologous variants between A, B and D, where A, B and D don't share the same allele","Inter-homoeologous variants between A, B and D bread wheat component genomes, where A, B and D don't share the same allele",1007);
INSERT INTO `variation_set` VALUES (8,"Inter-homoeologous variants between A, B and D, where A and B share the same allele","Inter-homoeologous variants between A, B and D bread wheat component genomes, where A and B share the same allele",1008);
INSERT INTO `variation_set` VALUES (9,"Inter-homoeologous variants between A, B and D, where A and D share the same allele","Inter-homoeologous variants between A, B and D bread wheat component genomes, where A and D share the same allele",1009);
INSERT INTO `variation_set` VALUES (10,"Inter-homoeologous variants between A, B and D, where B and D share the same allele","Inter-homoeologous variants between A, B and D bread wheat component genomes, where B and D share the same allele",1010);

INSERT INTO `variation_set` VALUES (11,"CerealsDB (all data)",1011);
INSERT INTO `variation_set` VALUES (12,"All variants classified as inter-homoeologues",1012);

-- variation_set_structure table

INSERT INTO `variation_set_structure` VALUES (11,1);
INSERT INTO `variation_set_structure` VALUES (11,2);
INSERT INTO `variation_set_structure` VALUES (11,3);
INSERT INTO `variation_set_structure` VALUES (12,4);
INSERT INTO `variation_set_structure` VALUES (12,5);
INSERT INTO `variation_set_structure` VALUES (12,6);
INSERT INTO `variation_set_structure` VALUES (12,7);
INSERT INTO `variation_set_structure` VALUES (12,8);
INSERT INTO `variation_set_structure` VALUES (12,9);
INSERT INTO `variation_set_structure` VALUES (12,10);
