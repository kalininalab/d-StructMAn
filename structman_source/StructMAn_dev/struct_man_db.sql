-- phpMyAdmin SQL Dump
-- version 4.0.10.11
-- http://www.phpmyadmin.net
--
-- Host: bioinfodb:3306
-- Generation Time: Aug 07, 2018 at 05:23 PM
-- Server version: 5.6.10
-- PHP Version: 5.6.17-0+deb8u1+mpi1

SET SQL_MODE = "NO_AUTO_VALUE_ON_ZERO";
SET time_zone = "+00:00";


/*!40101 SET @OLD_CHARACTER_SET_CLIENT=@@CHARACTER_SET_CLIENT */;
/*!40101 SET @OLD_CHARACTER_SET_RESULTS=@@CHARACTER_SET_RESULTS */;
/*!40101 SET @OLD_COLLATION_CONNECTION=@@COLLATION_CONNECTION */;
/*!40101 SET NAMES utf8 */;

--
-- Database: `struct_man_db_2`
--

-- --------------------------------------------------------

--
-- Table structure for table `Gene`
--

CREATE TABLE IF NOT EXISTS `Gene` (
  `Gene_Id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `Uniprot_Ac` varchar(64) DEFAULT NULL,
  `Genbank_Protein_Accession_Number` text,
  `Uniprot_Id` varchar(32) NOT NULL,
  `Original_Session` int(11) NOT NULL,
  `Error_Code` int(11) DEFAULT NULL,
  `Error` varchar(255) DEFAULT NULL,
  `Sequence` text,
  `Species` varchar(255) DEFAULT NULL,
  PRIMARY KEY (`Gene_Id`),
  UNIQUE KEY `Name` (`Uniprot_Ac`,`Uniprot_Id`),
  KEY `Uniprot_Ac` (`Uniprot_Ac`)
) ENGINE=InnoDB  DEFAULT CHARSET=latin1 AUTO_INCREMENT=101101 ;

-- --------------------------------------------------------

--
-- Table structure for table `GO_Term`
--

CREATE TABLE IF NOT EXISTS `GO_Term` (
  `GO_Term_Id` int(11) NOT NULL AUTO_INCREMENT,
  `Name` tinytext NOT NULL,
  `Id` varchar(64) NOT NULL,
  PRIMARY KEY (`GO_Term_Id`),
  UNIQUE KEY `Name` (`Id`)
) ENGINE=InnoDB  DEFAULT CHARSET=latin1 AUTO_INCREMENT=35390 ;

-- --------------------------------------------------------

--
-- Table structure for table `Ligand`
--

CREATE TABLE IF NOT EXISTS `Ligand` (
  `Ligand_Id` int(11) NOT NULL AUTO_INCREMENT,
  `Name` varchar(12) NOT NULL,
  `Smiles` text NOT NULL,
  `Inchi` text NOT NULL,
  PRIMARY KEY (`Ligand_Id`)
) ENGINE=InnoDB  DEFAULT CHARSET=latin1 AUTO_INCREMENT=15132 ;

-- --------------------------------------------------------

--
-- Table structure for table `Model`
--

CREATE TABLE IF NOT EXISTS `Model` (
  `Model_Id` int(11) NOT NULL AUTO_INCREMENT,
  `Name` varchar(64) NOT NULL,
  `DOPE_Score` float NOT NULL,
  `Mutation` int(11) NOT NULL,
  `Template` int(11) NOT NULL,
  PRIMARY KEY (`Model_Id`),
  UNIQUE KEY `Template` (`Template`),
  KEY `Mutation` (`Mutation`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1 AUTO_INCREMENT=1 ;

-- --------------------------------------------------------

--
-- Table structure for table `Mutation`
--

CREATE TABLE IF NOT EXISTS `Mutation` (
  `Mutation_Id` int(11) NOT NULL AUTO_INCREMENT,
  `Amino_Acid_Change` varchar(60) NOT NULL,
  `Gene` int(10) unsigned NOT NULL,
  `Polyphen2_Score` varchar(16) NOT NULL DEFAULT '-',
  PRIMARY KEY (`Mutation_Id`),
  KEY `Gene` (`Gene`)
) ENGINE=InnoDB  DEFAULT CHARSET=latin1 AUTO_INCREMENT=37328859 ;

-- --------------------------------------------------------

--
-- Table structure for table `Organism`
--

CREATE TABLE IF NOT EXISTS `Organism` (
  `Organism_Id` int(11) NOT NULL AUTO_INCREMENT,
  `Name` varchar(64) NOT NULL,
  PRIMARY KEY (`Organism_Id`),
  UNIQUE KEY `Name` (`Name`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1 AUTO_INCREMENT=1 ;

-- --------------------------------------------------------

--
-- Table structure for table `Pathway`
--

CREATE TABLE IF NOT EXISTS `Pathway` (
  `Pathway_Id` int(11) NOT NULL AUTO_INCREMENT,
  `Reactome_Id` varchar(32) NOT NULL,
  `Name` varchar(255) NOT NULL,
  PRIMARY KEY (`Pathway_Id`)
) ENGINE=InnoDB  DEFAULT CHARSET=latin1 AUTO_INCREMENT=3209 ;

-- --------------------------------------------------------

--
-- Table structure for table `RS_Annotation_Annotation`
--

CREATE TABLE IF NOT EXISTS `RS_Annotation_Annotation` (
  `Template_1` int(11) NOT NULL,
  `Template_2` int(11) NOT NULL,
  `Mutation_1` int(11) NOT NULL,
  `Mutation_2` int(11) NOT NULL,
  `Chain_1` char(1) NOT NULL,
  `Chain_2` char(1) NOT NULL,
  `Session` int(11) NOT NULL,
  `Distance` float NOT NULL,
  `Atompair` varchar(32) NOT NULL,
  KEY `Template_1` (`Template_1`),
  KEY `Template_2` (`Template_2`),
  KEY `Mutation_1` (`Mutation_1`),
  KEY `Mutation_2` (`Mutation_2`),
  KEY `Session` (`Session`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

-- --------------------------------------------------------

--
-- Table structure for table `RS_Annotation_Session`
--

CREATE TABLE IF NOT EXISTS `RS_Annotation_Session` (
  `Mutation` int(11) NOT NULL,
  `Template` int(11) NOT NULL,
  `Session` int(11) NOT NULL,
  KEY `Mutation` (`Mutation`),
  KEY `Template` (`Template`),
  KEY `Session` (`Session`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

-- --------------------------------------------------------

--
-- Table structure for table `RS_Gene_GO_Term`
--

CREATE TABLE IF NOT EXISTS `RS_Gene_GO_Term` (
  `Gene` int(10) unsigned NOT NULL,
  `GO_Term` int(11) NOT NULL,
  KEY `Gene` (`Gene`),
  KEY `GO_Term` (`GO_Term`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

-- --------------------------------------------------------

--
-- Table structure for table `RS_Gene_Pathway`
--

CREATE TABLE IF NOT EXISTS `RS_Gene_Pathway` (
  `Gene` int(10) unsigned NOT NULL,
  `Pathway` int(11) NOT NULL,
  KEY `Gene` (`Gene`),
  KEY `Pathway` (`Pathway`),
  KEY `Gene_2` (`Gene`),
  KEY `Gene_3` (`Gene`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

-- --------------------------------------------------------

--
-- Table structure for table `RS_Gene_Session`
--

CREATE TABLE IF NOT EXISTS `RS_Gene_Session` (
  `Gene` int(10) unsigned NOT NULL,
  `Session` int(11) NOT NULL,
  `Gene_Score` float DEFAULT NULL,
  KEY `Gene` (`Gene`,`Session`),
  KEY `Session` (`Session`),
  KEY `Gene_2` (`Gene`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

-- --------------------------------------------------------

--
-- Table structure for table `RS_Ligand_Template`
--

CREATE TABLE IF NOT EXISTS `RS_Ligand_Template` (
  `Ligand` int(11) NOT NULL,
  `Template` int(11) NOT NULL,
  `Chain` char(1) NOT NULL,
  `Residue` varchar(16) NOT NULL,
  KEY `Ligand` (`Ligand`,`Template`),
  KEY `Template` (`Template`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

-- --------------------------------------------------------

--
-- Table structure for table `RS_Mutation_Session`
--

CREATE TABLE IF NOT EXISTS `RS_Mutation_Session` (
  `Mutation` int(11) NOT NULL,
  `Session` int(11) NOT NULL,
  `New_AA` varchar(45) NOT NULL,
  `Tag` varchar(255) DEFAULT NULL,
  KEY `Mutation` (`Mutation`),
  KEY `Session` (`Session`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

-- --------------------------------------------------------

--
-- Table structure for table `RS_Mutation_Template`
--

CREATE TABLE IF NOT EXISTS `RS_Mutation_Template` (
  `Mutation` int(11) NOT NULL,
  `Template` int(11) NOT NULL,
  `Sub_Lig_Dist` text,
  `Sub_Chain_Distances` text,
  `Relative_Surface_Access` float DEFAULT NULL,
  `Secondary_Structure_Assignment` char(1) DEFAULT NULL,
  `Residue_Id` varchar(32) DEFAULT NULL,
  `Amino_Acid` char(1) DEFAULT NULL,
  `Candidate_Score` float DEFAULT NULL,
  `Homomer_Distances` varchar(1024) DEFAULT NULL,
  `Error_Code` int(11) DEFAULT NULL,
  `Error` varchar(255) DEFAULT NULL,
  `Ligand_Interaction_Degree` int(11) DEFAULT NULL,
  `Ligand_Interaction_Score` float DEFAULT NULL,
  `Chain_Interaction_Degree` int(11) DEFAULT NULL,
  `Chain_Interaction_Score` float DEFAULT NULL,
  `Short_Interaction_Degree` int(11) DEFAULT NULL,
  `Short_Interaction_Score` float DEFAULT NULL,
  `Medium_Interaction_Degree` int(11) DEFAULT NULL,
  `Medium_Interaction_Score` float DEFAULT NULL,
  `Long_Interaction_Degree` int(11) DEFAULT NULL,
  `Long_Interaction_Score` float DEFAULT NULL,
  KEY `Mutation` (`Mutation`),
  KEY `Template` (`Template`),
  KEY `Error` (`Error`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

-- --------------------------------------------------------

--
-- Table structure for table `Session`
--

CREATE TABLE IF NOT EXISTS `Session` (
  `Session_Id` int(11) NOT NULL AUTO_INCREMENT,
  `Input_File` varchar(255) NOT NULL,
  `Start` datetime NOT NULL,
  `End` datetime DEFAULT NULL,
  `Config_Mark` char(1) NOT NULL DEFAULT 'D',
  `Sequence_Identity_Threshold` float DEFAULT NULL,
  `Coverage_Threshold` float DEFAULT NULL,
  `Resolution_Threshold` float DEFAULT NULL,
  PRIMARY KEY (`Session_Id`)
) ENGINE=InnoDB  DEFAULT CHARSET=latin1 AUTO_INCREMENT=24 ;

-- --------------------------------------------------------

--
-- Table structure for table `Template`
--

CREATE TABLE IF NOT EXISTS `Template` (
  `Template_Id` int(11) NOT NULL AUTO_INCREMENT,
  `Name` varchar(64) NOT NULL,
  `Sequence_Identity` float NOT NULL,
  `Gene` int(11) unsigned NOT NULL,
  `Alignment_Length` float NOT NULL,
  `Target_Chain` char(1) NOT NULL,
  `Resolution` float NOT NULL,
  `R_Value` float NOT NULL,
  `Quality_Score` float NOT NULL,
  `Original_Target_Chain` char(1) NOT NULL,
  `Original_Chains` varchar(255) NOT NULL,
  `Chains` text NOT NULL,
  `Alignment` mediumtext NOT NULL,
  `Homooligomer` varchar(64) DEFAULT NULL,
  PRIMARY KEY (`Template_Id`),
  KEY `Gene` (`Gene`),
  KEY `Name` (`Name`),
  KEY `Template_Id` (`Template_Id`)
) ENGINE=InnoDB  DEFAULT CHARSET=latin1 AUTO_INCREMENT=2496892 ;

--
-- Constraints for dumped tables
--

--
-- Constraints for table `Model`
--
ALTER TABLE `Model`
  ADD CONSTRAINT `Model_ibfk_1` FOREIGN KEY (`Mutation`) REFERENCES `Mutation` (`Mutation_Id`) ON DELETE CASCADE ON UPDATE CASCADE,
  ADD CONSTRAINT `Model_ibfk_2` FOREIGN KEY (`Template`) REFERENCES `Template` (`Template_Id`) ON DELETE CASCADE ON UPDATE CASCADE;

--
-- Constraints for table `Mutation`
--
ALTER TABLE `Mutation`
  ADD CONSTRAINT `Mutation_ibfk_1` FOREIGN KEY (`Gene`) REFERENCES `Gene` (`Gene_Id`) ON DELETE CASCADE ON UPDATE CASCADE;

--
-- Constraints for table `RS_Annotation_Annotation`
--
ALTER TABLE `RS_Annotation_Annotation`
  ADD CONSTRAINT `RS_Annotation_Annotation_ibfk_1` FOREIGN KEY (`Template_1`) REFERENCES `Template` (`Template_Id`) ON DELETE CASCADE ON UPDATE CASCADE,
  ADD CONSTRAINT `RS_Annotation_Annotation_ibfk_2` FOREIGN KEY (`Template_2`) REFERENCES `Template` (`Template_Id`) ON DELETE CASCADE ON UPDATE CASCADE,
  ADD CONSTRAINT `RS_Annotation_Annotation_ibfk_3` FOREIGN KEY (`Mutation_1`) REFERENCES `Mutation` (`Mutation_Id`) ON DELETE CASCADE ON UPDATE CASCADE,
  ADD CONSTRAINT `RS_Annotation_Annotation_ibfk_4` FOREIGN KEY (`Mutation_2`) REFERENCES `Mutation` (`Mutation_Id`) ON DELETE CASCADE ON UPDATE CASCADE,
  ADD CONSTRAINT `RS_Annotation_Annotation_ibfk_5` FOREIGN KEY (`Session`) REFERENCES `Session` (`Session_Id`) ON DELETE CASCADE ON UPDATE CASCADE;

--
-- Constraints for table `RS_Annotation_Session`
--
ALTER TABLE `RS_Annotation_Session`
  ADD CONSTRAINT `RS_Annotation_Session_ibfk_1` FOREIGN KEY (`Mutation`) REFERENCES `Mutation` (`Mutation_Id`) ON DELETE CASCADE ON UPDATE CASCADE,
  ADD CONSTRAINT `RS_Annotation_Session_ibfk_2` FOREIGN KEY (`Template`) REFERENCES `Template` (`Template_Id`) ON DELETE CASCADE ON UPDATE CASCADE,
  ADD CONSTRAINT `RS_Annotation_Session_ibfk_3` FOREIGN KEY (`Session`) REFERENCES `Session` (`Session_Id`) ON DELETE CASCADE ON UPDATE CASCADE;

--
-- Constraints for table `RS_Gene_GO_Term`
--
ALTER TABLE `RS_Gene_GO_Term`
  ADD CONSTRAINT `RS_Gene_GO_Term_ibfk_1` FOREIGN KEY (`Gene`) REFERENCES `Gene` (`Gene_Id`) ON DELETE CASCADE ON UPDATE CASCADE,
  ADD CONSTRAINT `RS_Gene_GO_Term_ibfk_2` FOREIGN KEY (`GO_Term`) REFERENCES `GO_Term` (`GO_Term_Id`) ON DELETE CASCADE ON UPDATE CASCADE;

--
-- Constraints for table `RS_Gene_Pathway`
--
ALTER TABLE `RS_Gene_Pathway`
  ADD CONSTRAINT `RS_Gene_Pathway_ibfk_1` FOREIGN KEY (`Pathway`) REFERENCES `Pathway` (`Pathway_Id`) ON DELETE CASCADE ON UPDATE CASCADE,
  ADD CONSTRAINT `RS_Gene_Pathway_ibfk_2` FOREIGN KEY (`Gene`) REFERENCES `Gene` (`Gene_Id`) ON DELETE CASCADE ON UPDATE CASCADE;

--
-- Constraints for table `RS_Gene_Session`
--
ALTER TABLE `RS_Gene_Session`
  ADD CONSTRAINT `RS_Gene_Session_ibfk_1` FOREIGN KEY (`Session`) REFERENCES `Session` (`Session_Id`) ON DELETE CASCADE ON UPDATE CASCADE,
  ADD CONSTRAINT `RS_Gene_Session_ibfk_2` FOREIGN KEY (`Gene`) REFERENCES `Gene` (`Gene_Id`) ON DELETE CASCADE ON UPDATE CASCADE;

--
-- Constraints for table `RS_Ligand_Template`
--
ALTER TABLE `RS_Ligand_Template`
  ADD CONSTRAINT `RS_Ligand_Template_ibfk_1` FOREIGN KEY (`Ligand`) REFERENCES `Ligand` (`Ligand_Id`) ON DELETE CASCADE ON UPDATE CASCADE,
  ADD CONSTRAINT `RS_Ligand_Template_ibfk_2` FOREIGN KEY (`Template`) REFERENCES `Template` (`Template_Id`) ON DELETE CASCADE ON UPDATE CASCADE;

--
-- Constraints for table `RS_Mutation_Session`
--
ALTER TABLE `RS_Mutation_Session`
  ADD CONSTRAINT `RS_Mutation_Session_ibfk_1` FOREIGN KEY (`Mutation`) REFERENCES `Mutation` (`Mutation_Id`) ON DELETE CASCADE ON UPDATE CASCADE,
  ADD CONSTRAINT `RS_Mutation_Session_ibfk_2` FOREIGN KEY (`Session`) REFERENCES `Session` (`Session_Id`) ON DELETE CASCADE ON UPDATE CASCADE;

--
-- Constraints for table `RS_Mutation_Template`
--
ALTER TABLE `RS_Mutation_Template`
  ADD CONSTRAINT `RS_Mutation_Template_ibfk_1` FOREIGN KEY (`Mutation`) REFERENCES `Mutation` (`Mutation_Id`) ON DELETE CASCADE ON UPDATE CASCADE,
  ADD CONSTRAINT `RS_Mutation_Template_ibfk_2` FOREIGN KEY (`Template`) REFERENCES `Template` (`Template_Id`) ON DELETE CASCADE ON UPDATE CASCADE;

--
-- Constraints for table `Template`
--
ALTER TABLE `Template`
  ADD CONSTRAINT `Template_ibfk_1` FOREIGN KEY (`Gene`) REFERENCES `Gene` (`Gene_Id`) ON DELETE CASCADE ON UPDATE CASCADE;

/*!40101 SET CHARACTER_SET_CLIENT=@OLD_CHARACTER_SET_CLIENT */;
/*!40101 SET CHARACTER_SET_RESULTS=@OLD_CHARACTER_SET_RESULTS */;
/*!40101 SET COLLATION_CONNECTION=@OLD_COLLATION_CONNECTION */;
