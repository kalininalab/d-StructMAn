-- phpMyAdmin SQL Dump
-- version 4.9.1
-- https://www.phpmyadmin.net/
--
-- Host: wibi-guardian.helmholtz-hzi.de
-- Erstellungszeit: 17. Sep 2020 um 11:50
-- Server-Version: 10.1.46-MariaDB
-- PHP-Version: 7.0.33

SET SQL_MODE = "NO_AUTO_VALUE_ON_ZERO";
SET AUTOCOMMIT = 0;
START TRANSACTION;
SET time_zone = "+00:00";
SET @@auto_increment_increment=1;


/*!40101 SET @OLD_CHARACTER_SET_CLIENT=@@CHARACTER_SET_CLIENT */;
/*!40101 SET @OLD_CHARACTER_SET_RESULTS=@@CHARACTER_SET_RESULTS */;
/*!40101 SET @OLD_COLLATION_CONNECTION=@@COLLATION_CONNECTION */;
/*!40101 SET NAMES utf8mb4 */;

--
-- Database: `struct_man_db_1`
--
CREATE DATABASE IF NOT EXISTS `struct_man_db_1` DEFAULT CHARACTER SET latin1 COLLATE latin1_swedish_ci;
USE `struct_man_db_1`;
--

-- --------------------------------------------------------

--
-- Tabellenstruktur für Tabelle `Alignment`
--

CREATE TABLE `Alignment` (
  `Alignment_Id` int(11) UNSIGNED NOT NULL,
  `Gene` int(11) UNSIGNED NOT NULL,
  `Structure` int(11) UNSIGNED NOT NULL,
  `Sequence_Identity` float DEFAULT NULL,
  `Coverage` float DEFAULT NULL,
  `Alignment` mediumtext
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

-- --------------------------------------------------------

--
-- Tabellenstruktur für Tabelle `Complex`
--

CREATE TABLE `Complex` (
  `Complex_Id` int(11) UNSIGNED NOT NULL,
  `PDB` varchar(8) NOT NULL,
  `Resolution` float DEFAULT NULL,
  `Chains` text,
  `Homooligomers` text,
  `Ligand_Profile` text,
  `Metal_Profile` text,
  `Ion_Profile` text,
  `Chain_Chain_Profile` text
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

-- --------------------------------------------------------

--
-- Tabellenstruktur für Tabelle `Gene`
--

CREATE TABLE `Gene` (
  `Gene_Id` int(10) UNSIGNED NOT NULL,
  `Uniprot_Ac` varchar(255) DEFAULT NULL,
  `Genbank_Protein_Accession_Number` text DEFAULT NULL,
  `Uniprot_Id` varchar(32) DEFAULT NULL,
  `Original_Session` int(11) DEFAULT NULL,
  `Error_Code` int(11) DEFAULT NULL,
  `Error` varchar(255) DEFAULT NULL,
  `Sequence` text DEFAULT NULL,
  `Species` varchar(255) DEFAULT NULL
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

-- --------------------------------------------------------

--
-- Tabellenstruktur für Tabelle `GO_Term`
--

CREATE TABLE `GO_Term` (
  `GO_Term_Id` int(11) NOT NULL,
  `Name` tinytext NOT NULL,
  `Id` varchar(64) NOT NULL
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

-- --------------------------------------------------------

--
-- Tabellenstruktur für Tabelle `Indel`
--

CREATE TABLE `Indel` (
  `Indel_Id` int(10) UNSIGNED NOT NULL,
  `Wildtype_Protein` int(10) UNSIGNED NOT NULL,
  `Mutant_Protein` int(10) UNSIGNED NOT NULL,
  `Indel_Notation` text COLLATE utf8_unicode_ci NOT NULL
) ENGINE=InnoDB DEFAULT CHARSET=utf8 COLLATE=utf8_unicode_ci;

-- --------------------------------------------------------

--
-- Tabellenstruktur für Tabelle `Ligand`
--

CREATE TABLE `Ligand` (
  `Ligand_Id` int(11) NOT NULL,
  `Name` varchar(12) NOT NULL,
  `Smiles` text NOT NULL,
  `Inchi` text NOT NULL
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

-- --------------------------------------------------------

--
-- Tabellenstruktur für Tabelle `Mutation`
--

CREATE TABLE `Mutation` (
  `Mutation_Id` int(11) NOT NULL,
  `Amino_Acid_Change` varchar(60) DEFAULT NULL,
  `Residue_Id` varchar(8) DEFAULT NULL,
  `Gene` int(10) UNSIGNED DEFAULT NULL,
  `IUPRED` float DEFAULT NULL,
  `IUPRED_Glob` varchar(16) DEFAULT NULL,
  `Location` varchar(8) DEFAULT NULL,
  `Class` varchar(128) DEFAULT NULL,
  `RIN_Class` varchar(1024) DEFAULT NULL,
  `Simple_Class` varchar(64) DEFAULT NULL,
  `RIN_Simple_Class` varchar(64) DEFAULT NULL,
  `Interactions` text,
  `Confidence` float DEFAULT NULL,
  `Secondary_Structure` varchar(8) DEFAULT NULL,
  `Recommended_Structure` varchar(64) DEFAULT NULL,
  `Max_Seq_Structure` varchar(64) DEFAULT NULL,
  `Mapped_Structures` int(11) DEFAULT NULL,
  `RIN_Profile` text,
  `Modres_Score` float DEFAULT NULL,
  `Modres_Propensity` float DEFAULT NULL,
  `B_Factor` float DEFAULT NULL,
  `Weighted_Centrality_Scores` tinytext,
  `Weighted_Phi` float DEFAULT NULL,
  `Weighted_Psi` float DEFAULT NULL,
  `Intra_SSBOND_Propensity` float DEFAULT NULL,
  `Inter_SSBOND_Propensity` float DEFAULT NULL,
  `Intra_Link_Propensity` float DEFAULT NULL,
  `Inter_Link_Propensity` float DEFAULT NULL,
  `CIS_Conformation_Propensity` float DEFAULT NULL,
  `CIS_Follower_Propensity` float DEFAULT NULL,
  `Weighted_Inter_Chain_Median_KD` float DEFAULT NULL,
  `Weighted_Inter_Chain_Dist_Weighted_KD` float DEFAULT NULL,
  `Weighted_Inter_Chain_Median_RSA` float DEFAULT NULL,
  `Weighted_Inter_Chain_Dist_Weighted_RSA` float DEFAULT NULL,
  `Weighted_Intra_Chain_Median_KD` float DEFAULT NULL,
  `Weighted_Intra_Chain_Dist_Weighted_KD` float DEFAULT NULL,
  `Weighted_Intra_Chain_Median_RSA` float DEFAULT NULL,
  `Weighted_Intra_Chain_Dist_Weighted_RSA` float DEFAULT NULL
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

-- --------------------------------------------------------

--
-- Tabellenstruktur für Tabelle `Pathway`
--

CREATE TABLE `Pathway` (
  `Pathway_Id` int(11) NOT NULL,
  `Reactome_Id` varchar(32) NOT NULL,
  `Name` varchar(255) NOT NULL
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

-- --------------------------------------------------------

--
-- Tabellenstruktur für Tabelle `Residue`
--

CREATE TABLE `Residue` (
  `Residue_Id` int(11) UNSIGNED NOT NULL,
  `Structure` int(11) UNSIGNED NOT NULL,
  `Number` varchar(32) DEFAULT NULL,
  `Amino_Acid` char(1) DEFAULT NULL,
  `Sub_Lig_Dist` text,
  `Sub_Chain_Distances` text,
  `Relative_Surface_Access` float DEFAULT NULL,
  `Homomer_Distances` text,
  `Secondary_Structure_Assignment` char(1) DEFAULT NULL,
  `Interaction_Profile` varchar(1024) DEFAULT NULL,
  `Centralities` varchar(256) DEFAULT NULL,
  `B_Factor` float DEFAULT NULL,
  `Modres` tinyint(1) DEFAULT NULL,
  `PHI` float DEFAULT NULL,
  `PSI` float DEFAULT NULL,
  `Intra_SSBOND` tinyint(1) DEFAULT NULL,
  `SSBOND_Length` float DEFAULT NULL,
  `Intra_Link` tinyint(1) DEFAULT NULL,
  `Link_Length` float DEFAULT NULL,
  `CIS_Conformation` float DEFAULT NULL,
  `CIS_Follower` float DEFAULT NULL,
  `Inter_Chain_Median_KD` float DEFAULT NULL,
  `Inter_Chain_Dist_Weighted_KD` float DEFAULT NULL,
  `Inter_Chain_Median_RSA` float DEFAULT NULL,
  `Inter_Chain_Dist_Weighted_RSA` float DEFAULT NULL,
  `Intra_Chain_Median_KD` float DEFAULT NULL,
  `Intra_Chain_Dist_Weighted_KD` float DEFAULT NULL,
  `Intra_Chain_Median_RSA` float DEFAULT NULL,
  `Intra_Chain_Dist_Weighted_RSA` float DEFAULT NULL
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

-- --------------------------------------------------------

--
-- Tabellenstruktur für Tabelle `RS_Gene_GO_Term`
--

CREATE TABLE `RS_Gene_GO_Term` (
  `Gene` int(10) UNSIGNED NOT NULL,
  `GO_Term` int(11) NOT NULL
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

-- --------------------------------------------------------

--
-- Tabellenstruktur für Tabelle `RS_Gene_Pathway`
--

CREATE TABLE `RS_Gene_Pathway` (
  `Gene` int(10) UNSIGNED NOT NULL,
  `Pathway` int(11) NOT NULL
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

-- --------------------------------------------------------

--
-- Tabellenstruktur für Tabelle `RS_Gene_Session`
--

CREATE TABLE `RS_Gene_Session` (
  `Gene` int(10) UNSIGNED NOT NULL,
  `Session` int(11) NOT NULL,
  `Gene_Score` float DEFAULT NULL
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

-- --------------------------------------------------------

--
-- Tabellenstruktur für Tabelle `RS_Indel_Session`
--

CREATE TABLE `RS_Indel_Session` (
  `Session` int(11) NOT NULL,
  `Indel` int(10) UNSIGNED NOT NULL,
  `Tags` mediumtext COLLATE utf8_unicode_ci
) ENGINE=InnoDB DEFAULT CHARSET=utf8 COLLATE=utf8_unicode_ci;

-- --------------------------------------------------------

--
-- Tabellenstruktur für Tabelle `RS_Ligand_Structure`
--

CREATE TABLE `RS_Ligand_Structure` (
  `Ligand` int(11) NOT NULL,
  `Complex` int(11) UNSIGNED NOT NULL,
  `Chain` char(1) DEFAULT NULL,
  `Residue` varchar(16) DEFAULT NULL
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

-- --------------------------------------------------------

--
-- Tabellenstruktur für Tabelle `RS_Mutation_Session`
--

CREATE TABLE `RS_Mutation_Session` (
  `Mutation` int(11) NOT NULL,
  `Session` int(11) NOT NULL,
  `New_AA` varchar(45) NOT NULL,
  `Tag` mediumtext
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

-- --------------------------------------------------------

--
-- Tabellenstruktur für Tabelle `RS_Residue_Residue`
--

CREATE TABLE `RS_Residue_Residue` (
  `Residue_1` int(11) UNSIGNED NOT NULL,
  `Structure_1` int(11) UNSIGNED NOT NULL,
  `Residue_2` int(11) UNSIGNED NOT NULL,
  `Structure_2` int(11) UNSIGNED NOT NULL,
  `Distance` float NOT NULL
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

-- --------------------------------------------------------

--
-- Tabellenstruktur für Tabelle `Session`
--

CREATE TABLE `Session` (
  `Session_Id` int(11) NOT NULL,
  `Input_File` varchar(255) NOT NULL,
  `Start` datetime NOT NULL,
  `End` datetime DEFAULT NULL
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

-- --------------------------------------------------------

--
-- Tabellenstruktur für Tabelle `Structure`
--

CREATE TABLE `Structure` (
  `Structure_Id` int(11) UNSIGNED NOT NULL,
  `PDB` varchar(8) NOT NULL,
  `Chain` char(1) NOT NULL,
  `Homooligomer` varchar(256) DEFAULT NULL
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

--
-- Indizes der exportierten Tabellen
--

--
-- Indizes für die Tabelle `Alignment`
--
ALTER TABLE `Alignment`
  ADD PRIMARY KEY (`Alignment_Id`,`Gene`,`Structure`),
  ADD KEY `Gene` (`Gene`),
  ADD KEY `Structure` (`Structure`);

--
-- Indizes für die Tabelle `Complex`
--
ALTER TABLE `Complex`
  ADD PRIMARY KEY (`Complex_Id`);

--
-- Indizes für die Tabelle `Gene`
--
ALTER TABLE `Gene`
  ADD PRIMARY KEY (`Gene_Id`),
  ADD UNIQUE KEY `Name` (`Uniprot_Ac`,`Uniprot_Id`),
  ADD KEY `Uniprot_Ac` (`Uniprot_Ac`);

--
-- Indizes für die Tabelle `GO_Term`
--
ALTER TABLE `GO_Term`
  ADD PRIMARY KEY (`GO_Term_Id`),
  ADD UNIQUE KEY `Id` (`Id`) USING BTREE;

--
-- Indizes für die Tabelle `Indel`
--
ALTER TABLE `Indel`
  ADD PRIMARY KEY (`Indel_Id`),
  ADD KEY `Mutant_Protein` (`Mutant_Protein`),
  ADD KEY `Wildtype_Protein` (`Wildtype_Protein`);

--
-- Indizes für die Tabelle `Ligand`
--
ALTER TABLE `Ligand`
  ADD PRIMARY KEY (`Ligand_Id`);

--
-- Indizes für die Tabelle `Mutation`
--
ALTER TABLE `Mutation`
  ADD PRIMARY KEY (`Mutation_Id`),
  ADD KEY `Gene` (`Gene`);

--
-- Indizes für die Tabelle `Pathway`
--
ALTER TABLE `Pathway`
  ADD PRIMARY KEY (`Pathway_Id`);

--
-- Indizes für die Tabelle `Residue`
--
ALTER TABLE `Residue`
  ADD PRIMARY KEY (`Residue_Id`,`Structure`),
  ADD KEY `Structure` (`Structure`),
  ADD KEY `Residue_Id` (`Residue_Id`);

--
-- Indizes für die Tabelle `RS_Gene_GO_Term`
--
ALTER TABLE `RS_Gene_GO_Term`
  ADD KEY `Gene` (`Gene`),
  ADD KEY `GO_Term` (`GO_Term`);

--
-- Indizes für die Tabelle `RS_Gene_Pathway`
--
ALTER TABLE `RS_Gene_Pathway`
  ADD KEY `Gene` (`Gene`),
  ADD KEY `Pathway` (`Pathway`);

--
-- Indizes für die Tabelle `RS_Gene_Session`
--
ALTER TABLE `RS_Gene_Session`
  ADD KEY `Gene` (`Gene`,`Session`),
  ADD KEY `Session` (`Session`);

--
-- Indizes für die Tabelle `RS_Indel_Session`
--
ALTER TABLE `RS_Indel_Session`
  ADD KEY `Indel` (`Indel`),
  ADD KEY `Session` (`Session`);

--
-- Indizes für die Tabelle `RS_Ligand_Structure`
--
ALTER TABLE `RS_Ligand_Structure`
  ADD KEY `Ligand` (`Ligand`,`Complex`),
  ADD KEY `RS_Ligand_Structure_ibfk_2` (`Complex`);

--
-- Indizes für die Tabelle `RS_Mutation_Session`
--
ALTER TABLE `RS_Mutation_Session`
  ADD KEY `Mutation` (`Mutation`),
  ADD KEY `Session` (`Session`);

--
-- Indizes für die Tabelle `RS_Residue_Residue`
--
ALTER TABLE `RS_Residue_Residue`
  ADD KEY `Residue_1` (`Residue_1`),
  ADD KEY `Structure_1` (`Structure_1`),
  ADD KEY `Residue_2` (`Residue_2`),
  ADD KEY `Structure_2` (`Structure_2`);

--
-- Indizes für die Tabelle `Session`
--
ALTER TABLE `Session`
  ADD PRIMARY KEY (`Session_Id`);

--
-- Indizes für die Tabelle `Structure`
--
ALTER TABLE `Structure`
  ADD PRIMARY KEY (`Structure_Id`);

--
-- AUTO_INCREMENT für exportierte Tabellen
--

--
-- AUTO_INCREMENT für Tabelle `Alignment`
--
ALTER TABLE `Alignment`
  MODIFY `Alignment_Id` int(11) UNSIGNED NOT NULL AUTO_INCREMENT;

--
-- AUTO_INCREMENT für Tabelle `Complex`
--
ALTER TABLE `Complex`
  MODIFY `Complex_Id` int(11) UNSIGNED NOT NULL AUTO_INCREMENT;

--
-- AUTO_INCREMENT für Tabelle `Gene`
--
ALTER TABLE `Gene`
  MODIFY `Gene_Id` int(10) UNSIGNED NOT NULL AUTO_INCREMENT;

--
-- AUTO_INCREMENT für Tabelle `GO_Term`
--
ALTER TABLE `GO_Term`
  MODIFY `GO_Term_Id` int(11) NOT NULL AUTO_INCREMENT;

--
-- AUTO_INCREMENT für Tabelle `Indel`
--
ALTER TABLE `Indel`
  MODIFY `Indel_Id` int(10) UNSIGNED NOT NULL AUTO_INCREMENT;

--
-- AUTO_INCREMENT für Tabelle `Ligand`
--
ALTER TABLE `Ligand`
  MODIFY `Ligand_Id` int(11) NOT NULL AUTO_INCREMENT;

--
-- AUTO_INCREMENT für Tabelle `Mutation`
--
ALTER TABLE `Mutation`
  MODIFY `Mutation_Id` int(11) NOT NULL AUTO_INCREMENT;

--
-- AUTO_INCREMENT für Tabelle `Pathway`
--
ALTER TABLE `Pathway`
  MODIFY `Pathway_Id` int(11) NOT NULL AUTO_INCREMENT;

--
-- AUTO_INCREMENT für Tabelle `Residue`
--
ALTER TABLE `Residue`
  MODIFY `Residue_Id` int(11) UNSIGNED NOT NULL AUTO_INCREMENT;

--
-- AUTO_INCREMENT für Tabelle `Session`
--
ALTER TABLE `Session`
  MODIFY `Session_Id` int(11) NOT NULL AUTO_INCREMENT;

--
-- AUTO_INCREMENT für Tabelle `Structure`
--
ALTER TABLE `Structure`
  MODIFY `Structure_Id` int(11) UNSIGNED NOT NULL AUTO_INCREMENT;

--
-- Constraints der exportierten Tabellen
--

--
-- Constraints der Tabelle `Alignment`
--
ALTER TABLE `Alignment`
  ADD CONSTRAINT `Alignment_ibfk_1` FOREIGN KEY (`Gene`) REFERENCES `Gene` (`Gene_Id`) ON DELETE CASCADE ON UPDATE CASCADE,
  ADD CONSTRAINT `Alignment_ibfk_2` FOREIGN KEY (`Structure`) REFERENCES `Structure` (`Structure_Id`) ON DELETE CASCADE ON UPDATE CASCADE;

--
-- Constraints der Tabelle `Indel`
--
ALTER TABLE `Indel`
  ADD CONSTRAINT `Indel_ibfk_1` FOREIGN KEY (`Mutant_Protein`) REFERENCES `Gene` (`Gene_Id`) ON DELETE CASCADE ON UPDATE CASCADE,
  ADD CONSTRAINT `Indel_ibfk_2` FOREIGN KEY (`Wildtype_Protein`) REFERENCES `Gene` (`Gene_Id`) ON DELETE CASCADE ON UPDATE CASCADE;

--
-- Constraints der Tabelle `Mutation`
--
ALTER TABLE `Mutation`
  ADD CONSTRAINT `Mutation_ibfk_1` FOREIGN KEY (`Gene`) REFERENCES `Gene` (`Gene_Id`) ON DELETE CASCADE ON UPDATE CASCADE;

--
-- Constraints der Tabelle `Residue`
--
ALTER TABLE `Residue`
  ADD CONSTRAINT `Residue_ibfk_1` FOREIGN KEY (`Structure`) REFERENCES `Structure` (`Structure_Id`) ON DELETE CASCADE ON UPDATE CASCADE;

--
-- Constraints der Tabelle `RS_Gene_GO_Term`
--
ALTER TABLE `RS_Gene_GO_Term`
  ADD CONSTRAINT `RS_Gene_GO_Term_ibfk_1` FOREIGN KEY (`Gene`) REFERENCES `Gene` (`Gene_Id`) ON DELETE CASCADE ON UPDATE CASCADE,
  ADD CONSTRAINT `RS_Gene_GO_Term_ibfk_2` FOREIGN KEY (`GO_Term`) REFERENCES `GO_Term` (`GO_Term_Id`) ON DELETE CASCADE ON UPDATE CASCADE;

--
-- Constraints der Tabelle `RS_Gene_Pathway`
--
ALTER TABLE `RS_Gene_Pathway`
  ADD CONSTRAINT `RS_Gene_Pathway_ibfk_1` FOREIGN KEY (`Pathway`) REFERENCES `Pathway` (`Pathway_Id`) ON DELETE CASCADE ON UPDATE CASCADE,
  ADD CONSTRAINT `RS_Gene_Pathway_ibfk_2` FOREIGN KEY (`Gene`) REFERENCES `Gene` (`Gene_Id`) ON DELETE CASCADE ON UPDATE CASCADE;

--
-- Constraints der Tabelle `RS_Gene_Session`
--
ALTER TABLE `RS_Gene_Session`
  ADD CONSTRAINT `RS_Gene_Session_ibfk_1` FOREIGN KEY (`Session`) REFERENCES `Session` (`Session_Id`) ON DELETE CASCADE ON UPDATE CASCADE,
  ADD CONSTRAINT `RS_Gene_Session_ibfk_2` FOREIGN KEY (`Gene`) REFERENCES `Gene` (`Gene_Id`) ON DELETE CASCADE ON UPDATE CASCADE;

--
-- Constraints der Tabelle `RS_Indel_Session`
--
ALTER TABLE `RS_Indel_Session`
  ADD CONSTRAINT `RS_Indel_Session_ibfk_1` FOREIGN KEY (`Indel`) REFERENCES `Indel` (`Indel_Id`) ON DELETE CASCADE ON UPDATE CASCADE,
  ADD CONSTRAINT `RS_Indel_Session_ibfk_2` FOREIGN KEY (`Session`) REFERENCES `Session` (`Session_Id`) ON DELETE CASCADE ON UPDATE CASCADE;

--
-- Constraints der Tabelle `RS_Ligand_Structure`
--
ALTER TABLE `RS_Ligand_Structure`
  ADD CONSTRAINT `RS_Ligand_Structure_ibfk_1` FOREIGN KEY (`Ligand`) REFERENCES `Ligand` (`Ligand_Id`) ON DELETE CASCADE ON UPDATE CASCADE,
  ADD CONSTRAINT `RS_Ligand_Structure_ibfk_2` FOREIGN KEY (`Complex`) REFERENCES `Complex` (`Complex_Id`) ON DELETE CASCADE ON UPDATE CASCADE;

--
-- Constraints der Tabelle `RS_Mutation_Session`
--
ALTER TABLE `RS_Mutation_Session`
  ADD CONSTRAINT `RS_Mutation_Session_ibfk_1` FOREIGN KEY (`Mutation`) REFERENCES `Mutation` (`Mutation_Id`) ON DELETE CASCADE ON UPDATE CASCADE,
  ADD CONSTRAINT `RS_Mutation_Session_ibfk_2` FOREIGN KEY (`Session`) REFERENCES `Session` (`Session_Id`) ON DELETE CASCADE ON UPDATE CASCADE;

--
-- Constraints der Tabelle `RS_Residue_Residue`
--
ALTER TABLE `RS_Residue_Residue`
  ADD CONSTRAINT `RS_Residue_Residue_ibfk_1` FOREIGN KEY (`Residue_1`) REFERENCES `Residue` (`Residue_Id`) ON DELETE CASCADE ON UPDATE CASCADE,
  ADD CONSTRAINT `RS_Residue_Residue_ibfk_2` FOREIGN KEY (`Structure_1`) REFERENCES `Structure` (`Structure_Id`) ON DELETE CASCADE ON UPDATE CASCADE,
  ADD CONSTRAINT `RS_Residue_Residue_ibfk_3` FOREIGN KEY (`Residue_2`) REFERENCES `Residue` (`Residue_Id`) ON DELETE CASCADE ON UPDATE CASCADE,
  ADD CONSTRAINT `RS_Residue_Residue_ibfk_4` FOREIGN KEY (`Structure_2`) REFERENCES `Structure` (`Structure_Id`) ON DELETE CASCADE ON UPDATE CASCADE;
COMMIT;

/*!40101 SET CHARACTER_SET_CLIENT=@OLD_CHARACTER_SET_CLIENT */;
/*!40101 SET CHARACTER_SET_RESULTS=@OLD_CHARACTER_SET_RESULTS */;
/*!40101 SET COLLATION_CONNECTION=@OLD_COLLATION_CONNECTION */;
