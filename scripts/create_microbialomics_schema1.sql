-- MySQL dump 10.13  Distrib 5.6.13, for osx10.7 (x86_64)
--
-- Host: localhost    Database: microbialomics_v3
-- ------------------------------------------------------
-- Server version	5.6.13

/*!40101 SET @OLD_CHARACTER_SET_CLIENT=@@CHARACTER_SET_CLIENT */;
/*!40101 SET @OLD_CHARACTER_SET_RESULTS=@@CHARACTER_SET_RESULTS */;
/*!40101 SET @OLD_COLLATION_CONNECTION=@@COLLATION_CONNECTION */;
/*!40101 SET NAMES utf8 */;
/*!40103 SET @OLD_TIME_ZONE=@@TIME_ZONE */;
/*!40103 SET TIME_ZONE='+00:00' */;
/*!40014 SET @OLD_UNIQUE_CHECKS=@@UNIQUE_CHECKS, UNIQUE_CHECKS=0 */;
/*!40014 SET @OLD_FOREIGN_KEY_CHECKS=@@FOREIGN_KEY_CHECKS, FOREIGN_KEY_CHECKS=0 */;
/*!40101 SET @OLD_SQL_MODE=@@SQL_MODE, SQL_MODE='NO_AUTO_VALUE_ON_ZERO' */;
/*!40111 SET @OLD_SQL_NOTES=@@SQL_NOTES, SQL_NOTES=0 */;

--
-- Table structure for table `GENOME`
--

DROP TABLE IF EXISTS `GENOME`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `GENOME` (
  `genome_id` smallint(5) unsigned NOT NULL AUTO_INCREMENT,
  `genome_constant_id` varchar(100) NOT NULL,
  `genome_name` varchar(100) NOT NULL,
  `genome_description` varchar(255) DEFAULT NULL,
  `genus_species` varchar(100) NOT NULL,
  `genome_length` int(10) unsigned NOT NULL,
  `genome_sequence` longtext NOT NULL,
  `sequencing_date` date DEFAULT NULL,
  `pmid` int(10) unsigned ,
  `sequencing_method` varchar(255) ,
  `assembly_method` varchar(255) ,
  `aerobic_tolerance` varchar(255) ,
  `growth_conditions` varchar(255) ,
  `public_genome` tinyint(1) ,
  PRIMARY KEY (`genome_id`),
  KEY `genome_name` (`genome_name`),
  KEY `genus_species` (`genus_species`),
  KEY `genome_id` (`genome_id`)
) ENGINE=MyISAM  DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `build`
--

DROP TABLE IF EXISTS `build`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `build` (
  `build_id` int(11) NOT NULL,
  `build_name` varchar(32) NOT NULL,
  `build_date` date NOT NULL,
  `description` text NOT NULL
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `contig`
--

DROP TABLE IF EXISTS `contig`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `contig` (
  `contig_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `genome_id` smallint(5) unsigned NOT NULL,
  `contig_description` varchar(255) DEFAULT NULL,
  `start_pos` int(11) NOT NULL,
  `stop_pos` int(11) NOT NULL,
  PRIMARY KEY (`contig_id`),
  KEY `genome_id` (`genome_id`)
) ENGINE=MyISAM  DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `feature`
--

DROP TABLE IF EXISTS `feature`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `feature` (
  `feature_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `genome_id` smallint(5) unsigned NOT NULL,
  `genbank_id` int unsigned ,
  `feature_name` varchar(64) NOT NULL,
  `feature_symbol` varchar(64) NOT NULL,
  `feature_description` varchar(255) ,
  `feature_type` enum('CDS','TU','5S rRNA','16S rRNA','23S rRNA','miRNA','transcription start site','operon','tRNA','sRNA','miscRNA','TransTerm') NOT NULL,
  `feature_evidence` varchar(64) ,
  `feature_score` float ,
  `start_pos` int(11) NOT NULL,
  `stop_pos` int(11) NOT NULL,
  `direction` enum('F','R') NOT NULL,
  `primary_feature` enum('Y','N') NOT NULL,
  PRIMARY KEY (`feature_id`),
  KEY `feature_name` (`feature_name`),
  KEY `feature_symbol` (`feature_symbol`),
  KEY `feature_type` (`feature_type`),
  KEY `primary_feature` (`primary_feature`),
  KEY `feature_id` (`feature_id`,`genome_id`),
  KEY `genome_id` (`genome_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `feature_to_function`
--

DROP TABLE IF EXISTS `feature_to_function`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `feature_to_function` (
  `feature_id` int(10) unsigned NOT NULL,
  `function_id` int(10) unsigned NOT NULL,
  `score` float DEFAULT NULL,
  `eval` double DEFAULT NULL,
  `pval` double DEFAULT NULL,
  `direct_annotation` enum('Y','N') NOT NULL,
  `best_hit` enum('Y','N') DEFAULT 'N',
  UNIQUE KEY `feature_id_2` (`feature_id`,`function_id`),
  KEY `feature_id` (`feature_id`),
  KEY `function_id` (`function_id`),
  KEY `best_hit` (`best_hit`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `function`
--

DROP TABLE IF EXISTS `function`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `function` (
  `function_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `function_accession` varchar(100) NOT NULL,
  `function_name` varchar(255) NOT NULL,
  `function_description` varchar(1000) DEFAULT NULL,
  `function_type` enum('COG category','COG (STRING v7.1)','human_defined','expression_modulation','TIGRFAM','PFAM','GO','CELLO','CAZy','KEGG','KEGG category 1','KEGG category 2','KEGG pathway') NOT NULL,
  PRIMARY KEY (`function_id`)
) ENGINE=MyISAM  DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `genetics`
--

DROP TABLE IF EXISTS `genetics`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `genetics` (
  `feature_id` int(10) unsigned NOT NULL,
  `perturbation_type` enum('knockout','overexpression') NOT NULL,
  `hit_type` enum('direct','polar','NA') DEFAULT NULL,
  `strain_id` varchar(256) NOT NULL,
  `description` varchar(256) DEFAULT NULL,
  `multiple_matches` enum('Y','N') DEFAULT NULL,
  `insertion_pos` int(11) DEFAULT NULL,
  KEY `strain_id` (`strain_id`),
  KEY `feature_id` (`feature_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `model_species_homologs`
--

DROP TABLE IF EXISTS `model_species_homologs`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `model_species_homologs` (
  `feature_id_model` int(10) unsigned NOT NULL,
  `feature_id_homolog` int(10) unsigned NOT NULL,
  KEY `feature_id_homolog` (`feature_id_homolog`),
  KEY `feature_id_model` (`feature_id_model`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `protein_seq`
--

DROP TABLE IF EXISTS `protein_seq`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `protein_seq` (
  `feature_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `protein_sequence` text NOT NULL,
  PRIMARY KEY (`feature_id`)
) ENGINE=MyISAM  DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `regulation`
--

DROP TABLE IF EXISTS `regulation`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `regulation` (
  `feature_id_regulator` int(10) unsigned NOT NULL,
  `feature_id_target` int(10) unsigned NOT NULL,
  `regulation_type` enum('protein-protein','transcription factor-promoter') NOT NULL,
  `evidence` enum('ChIP','expression','yeast','footprint') NOT NULL,
  `source` enum('RegulonDB') NOT NULL,
  KEY `feature_id_regulator` (`feature_id_regulator`),
  KEY `feature_id_target` (`feature_id_target`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;
/*!40103 SET TIME_ZONE=@OLD_TIME_ZONE */;

/*!40101 SET SQL_MODE=@OLD_SQL_MODE */;
/*!40014 SET FOREIGN_KEY_CHECKS=@OLD_FOREIGN_KEY_CHECKS */;
/*!40014 SET UNIQUE_CHECKS=@OLD_UNIQUE_CHECKS */;
/*!40101 SET CHARACTER_SET_CLIENT=@OLD_CHARACTER_SET_CLIENT */;
/*!40101 SET CHARACTER_SET_RESULTS=@OLD_CHARACTER_SET_RESULTS */;
/*!40101 SET COLLATION_CONNECTION=@OLD_COLLATION_CONNECTION */;
/*!40111 SET SQL_NOTES=@OLD_SQL_NOTES */;

-- Dump completed on 2014-07-10 11:49:25
