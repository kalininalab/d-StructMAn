-- phpMyAdmin SQL Dump
-- version 4.0.10.11
-- http://www.phpmyadmin.net
--
-- Host: bioinfodb:3306
-- Generation Time: Feb 24, 2020 at 05:11 PM
-- Server version: 5.6.10
-- PHP Version: 5.6.33-0+deb8u1+mpi1

SET SQL_MODE = "NO_AUTO_VALUE_ON_ZERO";
SET time_zone = "+00:00";


/*!40101 SET @OLD_CHARACTER_SET_CLIENT=@@CHARACTER_SET_CLIENT */;
/*!40101 SET @OLD_CHARACTER_SET_RESULTS=@@CHARACTER_SET_RESULTS */;
/*!40101 SET @OLD_COLLATION_CONNECTION=@@COLLATION_CONNECTION */;
/*!40101 SET NAMES utf8 */;

--
-- Database: `struct_man_db_uniprot`
--

-- --------------------------------------------------------

--
-- Table structure for table `AC_ID`
--

CREATE TABLE IF NOT EXISTS `AC_ID` (
  `AC_ID_Id` int(11) NOT NULL AUTO_INCREMENT,
  `Uniprot_Ac` varchar(16) NOT NULL,
  `Ac_Hash` varchar(2) NOT NULL,
  `Uniprot_Id` varchar(32) NOT NULL,
  `Id_Hash` varchar(2) NOT NULL,
  PRIMARY KEY (`AC_ID_Id`),
  KEY `Uniprot_Ac` (`Uniprot_Ac`,`Ac_Hash`,`Uniprot_Id`,`Id_Hash`),
  KEY `Ac_Hash` (`Ac_Hash`),
  KEY `Uniprot_Id` (`Uniprot_Id`),
  KEY `Id_Hash` (`Id_Hash`)
) ENGINE=InnoDB  DEFAULT CHARSET=latin1;

-- --------------------------------------------------------

--
-- Table structure for table `AC_Refseq`
--

CREATE TABLE IF NOT EXISTS `AC_Refseq` (
  `AC_Refseq_Id` int(11) NOT NULL AUTO_INCREMENT,
  `Uniprot_Ac` varchar(16) NOT NULL,
  `Ac_Hash` varchar(2) NOT NULL,
  `Refseq` varchar(16) NOT NULL,
  `Refseq_Hash` varchar(2) NOT NULL,
  PRIMARY KEY (`AC_Refseq_Id`),
  KEY `Uniprot_Ac` (`Uniprot_Ac`,`Ac_Hash`,`Refseq`,`Refseq_Hash`),
  KEY `Refseq` (`Refseq`),
  KEY `Refseq_Hash` (`Refseq_Hash`),
  KEY `Ac_Hash` (`Ac_Hash`)
) ENGINE=InnoDB  DEFAULT CHARSET=latin1;

-- --------------------------------------------------------

--
-- Table structure for table `Sequences`
--

CREATE TABLE IF NOT EXISTS `Sequences` (
  `Sequence_Id` int(11) NOT NULL AUTO_INCREMENT,
  `Uniprot_Ac` varchar(32) NOT NULL,
  `Sequence` text,
  `Disorder_Scores` text,
  `Disorder_Regions` text,
  PRIMARY KEY (`Sequence_Id`,`Uniprot_Ac`),
  KEY `Uniprot_Ac` (`Uniprot_Ac`)
) ENGINE=InnoDB  DEFAULT CHARSET=latin1;

/*!40101 SET CHARACTER_SET_CLIENT=@OLD_CHARACTER_SET_CLIENT */;
/*!40101 SET CHARACTER_SET_RESULTS=@OLD_CHARACTER_SET_RESULTS */;
/*!40101 SET COLLATION_CONNECTION=@OLD_COLLATION_CONNECTION */;
