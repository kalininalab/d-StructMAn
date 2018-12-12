-- phpMyAdmin SQL Dump
-- version 4.0.10.11
-- http://www.phpmyadmin.net
--
-- Host: bioinfodb:3306
-- Generation Time: Dec 12, 2018 at 12:54 PM
-- Server version: 5.6.10
-- PHP Version: 5.6.33-0+deb8u1+mpi1

SET SQL_MODE = "NO_AUTO_VALUE_ON_ZERO";
SET time_zone = "+00:00";


/*!40101 SET @OLD_CHARACTER_SET_CLIENT=@@CHARACTER_SET_CLIENT */;
/*!40101 SET @OLD_CHARACTER_SET_RESULTS=@@CHARACTER_SET_RESULTS */;
/*!40101 SET @OLD_COLLATION_CONNECTION=@@COLLATION_CONNECTION */;
/*!40101 SET NAMES utf8 */;

--
-- Database: `struct_man_db_1`
--
CREATE DATABASE IF NOT EXISTS `struct_man_db_1` DEFAULT CHARACTER SET latin1 COLLATE latin1_swedish_ci;
USE `struct_man_db_1`;

-- --------------------------------------------------------

--
-- Table structure for table `Alignment`
--

CREATE TABLE IF NOT EXISTS `Alignment` (
  `Alignment_Id` int(11) unsigned NOT NULL AUTO_INCREMENT,
  `Gene` int(11) unsigned NOT NULL,
  `Structure` int(11) unsigned NOT NULL,
  `Sequence_Identity` float NOT NULL,
  `Coverage` float NOT NULL,
  `Alignment` mediumtext,
  PRIMARY KEY (`Alignment_Id`,`Gene`,`Structure`),
  KEY `Gene` (`Gene`),
  KEY `Structure` (`Structure`)
) ENGINE=InnoDB  DEFAULT CHARSET=latin1 AUTO_INCREMENT=1377734 ;

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
) ENGINE=InnoDB  DEFAULT CHARSET=latin1 AUTO_INCREMENT=48210 ;

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
) ENGINE=InnoDB  DEFAULT CHARSET=latin1 AUTO_INCREMENT=23597 ;

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
) ENGINE=InnoDB  DEFAULT CHARSET=latin1 AUTO_INCREMENT=15883 ;

-- --------------------------------------------------------

--
-- Table structure for table `Mutation`
--

CREATE TABLE IF NOT EXISTS `Mutation` (
  `Mutation_Id` int(11) NOT NULL AUTO_INCREMENT,
  `Amino_Acid_Change` varchar(60) NOT NULL,
  `Gene` int(10) unsigned NOT NULL,
  `IUPRED` float DEFAULT NULL,
  `IUPRED_Glob` tinyint(1) DEFAULT NULL,
  PRIMARY KEY (`Mutation_Id`),
  KEY `Gene` (`Gene`)
) ENGINE=InnoDB  DEFAULT CHARSET=latin1 AUTO_INCREMENT=19419541 ;

-- --------------------------------------------------------

--
-- Table structure for table `Pathway`
--

CREATE TABLE IF NOT EXISTS `Pathway` (
  `Pathway_Id` int(11) NOT NULL AUTO_INCREMENT,
  `Reactome_Id` varchar(32) NOT NULL,
  `Name` varchar(255) NOT NULL,
  PRIMARY KEY (`Pathway_Id`)
) ENGINE=InnoDB  DEFAULT CHARSET=latin1 AUTO_INCREMENT=2351 ;

-- --------------------------------------------------------

--
-- Table structure for table `Residue`
--

CREATE TABLE IF NOT EXISTS `Residue` (
  `Residue_Id` int(11) unsigned NOT NULL AUTO_INCREMENT,
  `Structure` int(11) unsigned NOT NULL,
  `Number` varchar(32) DEFAULT NULL,
  `Amino_Acid` char(1) DEFAULT NULL,
  `Sub_Lig_Dist` text,
  `Sub_Chain_Distances` text,
  `Relative_Surface_Access` float DEFAULT NULL,
  `Homomer_Distances` text,
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
  `Secondary_Structure_Assignment` char(1) DEFAULT NULL,
  PRIMARY KEY (`Residue_Id`,`Structure`),
  KEY `Structure` (`Structure`),
  KEY `Residue_Id` (`Residue_Id`)
) ENGINE=InnoDB  DEFAULT CHARSET=latin1 AUTO_INCREMENT=25064386 ;

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
-- Table structure for table `RS_Ligand_Structure`
--

CREATE TABLE IF NOT EXISTS `RS_Ligand_Structure` (
  `Ligand` int(11) NOT NULL,
  `Structure` int(11) unsigned NOT NULL,
  `Chain` char(1) DEFAULT NULL,
  `Residue` varchar(16) DEFAULT NULL,
  PRIMARY KEY (`Ligand`,`Structure`),
  KEY `Structure` (`Structure`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

-- --------------------------------------------------------

--
-- Table structure for table `RS_Mutation_Residue`
--

CREATE TABLE IF NOT EXISTS `RS_Mutation_Residue` (
  `Mutation` int(11) NOT NULL,
  `Residue` int(11) unsigned NOT NULL,
  `Structure` int(11) unsigned NOT NULL,
  `Gene` int(10) unsigned NOT NULL,
  KEY `Structure` (`Structure`,`Gene`),
  KEY `Gene` (`Gene`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

-- --------------------------------------------------------

--
-- Table structure for table `RS_Mutation_Session`
--

CREATE TABLE IF NOT EXISTS `RS_Mutation_Session` (
  `Mutation` int(11) NOT NULL,
  `Session` int(11) NOT NULL,
  `New_AA` varchar(45) NOT NULL,
  `Tag` text,
  KEY `Mutation` (`Mutation`),
  KEY `Session` (`Session`)
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
  PRIMARY KEY (`Session_Id`)
) ENGINE=InnoDB  DEFAULT CHARSET=latin1 AUTO_INCREMENT=3 ;

-- --------------------------------------------------------

--
-- Table structure for table `Structure`
--

CREATE TABLE IF NOT EXISTS `Structure` (
  `Structure_Id` int(11) unsigned NOT NULL AUTO_INCREMENT,
  `PDB` varchar(8) NOT NULL,
  `Chain` char(1) NOT NULL,
  `Resolution` float NOT NULL,
  `Homooligomer` varchar(64) DEFAULT NULL,
  `Chains` text,
  PRIMARY KEY (`Structure_Id`)
) ENGINE=InnoDB  DEFAULT CHARSET=latin1 AUTO_INCREMENT=95427 ;

--
-- Constraints for dumped tables
--

--
-- Constraints for table `Alignment`
--
ALTER TABLE `Alignment`
  ADD CONSTRAINT `Alignment_ibfk_1` FOREIGN KEY (`Gene`) REFERENCES `Gene` (`Gene_Id`) ON DELETE CASCADE ON UPDATE CASCADE,
  ADD CONSTRAINT `Alignment_ibfk_2` FOREIGN KEY (`Structure`) REFERENCES `Structure` (`Structure_Id`) ON DELETE CASCADE ON UPDATE CASCADE;

--
-- Constraints for table `Mutation`
--
ALTER TABLE `Mutation`
  ADD CONSTRAINT `Mutation_ibfk_1` FOREIGN KEY (`Gene`) REFERENCES `Gene` (`Gene_Id`) ON DELETE CASCADE ON UPDATE CASCADE;

--
-- Constraints for table `Residue`
--
ALTER TABLE `Residue`
  ADD CONSTRAINT `Residue_ibfk_1` FOREIGN KEY (`Structure`) REFERENCES `Structure` (`Structure_Id`) ON DELETE CASCADE ON UPDATE CASCADE;

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
-- Constraints for table `RS_Ligand_Structure`
--
ALTER TABLE `RS_Ligand_Structure`
  ADD CONSTRAINT `RS_Ligand_Structure_ibfk_1` FOREIGN KEY (`Structure`) REFERENCES `Structure` (`Structure_Id`) ON DELETE CASCADE ON UPDATE CASCADE,
  ADD CONSTRAINT `RS_Ligand_Structure_ibfk_2` FOREIGN KEY (`Ligand`) REFERENCES `Ligand` (`Ligand_Id`) ON DELETE CASCADE ON UPDATE CASCADE;

--
-- Constraints for table `RS_Mutation_Residue`
--
ALTER TABLE `RS_Mutation_Residue`
  ADD CONSTRAINT `RS_Mutation_Residue_ibfk_3` FOREIGN KEY (`Structure`) REFERENCES `Structure` (`Structure_Id`) ON DELETE CASCADE ON UPDATE CASCADE,
  ADD CONSTRAINT `RS_Mutation_Residue_ibfk_4` FOREIGN KEY (`Gene`) REFERENCES `Gene` (`Gene_Id`) ON DELETE CASCADE ON UPDATE CASCADE;

--
-- Constraints for table `RS_Mutation_Session`
--
ALTER TABLE `RS_Mutation_Session`
  ADD CONSTRAINT `RS_Mutation_Session_ibfk_1` FOREIGN KEY (`Mutation`) REFERENCES `Mutation` (`Mutation_Id`) ON DELETE CASCADE ON UPDATE CASCADE,
  ADD CONSTRAINT `RS_Mutation_Session_ibfk_2` FOREIGN KEY (`Session`) REFERENCES `Session` (`Session_Id`) ON DELETE CASCADE ON UPDATE CASCADE;

/*!40101 SET CHARACTER_SET_CLIENT=@OLD_CHARACTER_SET_CLIENT */;
/*!40101 SET CHARACTER_SET_RESULTS=@OLD_CHARACTER_SET_RESULTS */;
/*!40101 SET COLLATION_CONNECTION=@OLD_COLLATION_CONNECTION */;
