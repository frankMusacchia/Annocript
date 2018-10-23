/*
*   dbCreation.sql
*   author: Francesco Musacchia
*   summary: this script is needed forse parse_swissprot
*            script and creates all the tables and a view
* ################################
* BE CAREFUL!!!!!
* When you update this script to change tables in the database,
* you must also change the version here and the variable $dbTablesVersion
* in annocript.pl
* #################################
*/


/*FOREIGN KEYS HAVE TO BE WRITTEN EVERYTIME AT THE END OF
THE TABLE DEFINITION
PRAGMA foreign_keys = ON;/*In SQLITE foreign_keys are de-activated by default*/

/*This table contains indication of the database where the proteins have been taken (Swiss-Prot or Trembl)*/
CREATE TABLE dbType (
    id TINYINT unsigned PRIMARY KEY, 
    type varchar(20) not null,
    unique (type) );

/*Descriptions of the proteins*/    
CREATE TABLE uniprotkbDescType (
    id MEDIUMINT unsigned PRIMARY KEY,
    type MEDIUMTEXT not null);

/*Organisms*/
CREATE TABLE OSType (
    id MEDIUMINT unsigned PRIMARY KEY,
    type TEXT not null);
    
CREATE TABLE GNType (
    id MEDIUMINT unsigned PRIMARY KEY,
    type MEDIUMTEXT  not null);
    
CREATE TABLE PEType (
    id TINYINT unsigned PRIMARY KEY,
    type  varchar(10) not null,
    unique (type));
    
CREATE TABLE SVType (
    id TINYINT unsigned PRIMARY KEY,
    type  varchar(10) not null,
    unique (type));

/*This table will contain the GO ids with their descriptions*/    
CREATE TABLE goTable (
    id INT unsigned PRIMARY KEY,
    go_id  varchar(50) not null,
    definition MEDIUMTEXT,
    division ENUM('C','F','P'),
    unique (go_id) );
 
 /*Domains with descriptions*/
CREATE TABLE cdTable (
    cd_id INT unsigned PRIMARY KEY,
    genericId  varchar(50) not null,
    cdName  varchar(50) not null,
    cdShortName varchar(50) not null,
    cdDescription MEDIUMTEXT,
    pssmLength INT unsigned,
    unique (cdName) );

/*A table with the mapping from the domains to their GO ids*/   
CREATE TABLE cdMappingsTable (
    cd_id INT unsigned not null,
    go_ids MEDIUMTEXT,
    FOREIGN KEY (cd_id) REFERENCES cdTable(cd_id),
    PRIMARY KEY (cd_id)
    );    
   
/*This table will contain the EC Enzyme ids with their descriptions*/        
CREATE TABLE enzymeTable (
	id INT unsigned PRIMARY KEY, 
	enz_id varchar(50) not null, 
	description mediumtext, 
	unique (enz_id));

/*This table will contain the Pathways ids with their descriptions*/        
CREATE TABLE pathwaysTable (
	id INT unsigned PRIMARY KEY, 
	description mediumtext, 
	unique (id));

/*This is a view that permits to access the GO terms ids of a 
given doomain using its ID*/   
create view cdViewGO as
    select c.cdName as cdName, m.go_ids as go_ids
        FROM  cdMappingsTable m LEFT JOIN cdTable c ON  m.cd_id = c.cd_id;
          
/*This is the main table with the proteins and all the ids referred to descriptions in other
tables*/  
CREATE TABLE uniprotkb (
            id INT unsigned PRIMARY KEY,
            trId varchar (100) not null,
            seqId varchar (50) not null,
            dbId TINYINT unsigned not null,
            descId MEDIUMINT unsigned not null,
            OSId MEDIUMINT  unsigned not null,
            GNId MEDIUMINT unsigned not null,
            PEId TINYINT  unsigned not null,
            SVId TINYINT unsigned not null,
            FOREIGN KEY(dbId) REFERENCES dbType(id),
            FOREIGN KEY(descId) REFERENCES uniprotkbDescType(id),
            FOREIGN KEY(OSId) REFERENCES OSType(id),
            FOREIGN KEY(GNId) REFERENCES GNType(id),
            FOREIGN KEY(PEId) REFERENCES PEType(id),
            FOREIGN KEY(SVId) REFERENCES SVType(id),
            /*FOREIGN KEY(id) REFERENCES mappingsTable(uniprot_id),*/
            unique(seqId)
            );
                                        
/*A table with the mapping from the proteins to their GO and ENZYME ids*/   
CREATE TABLE uniprotkbMappingsTable (
    uniprot_id INT unsigned not null,
    go_ids MEDIUMTEXT,
    enzyme_ids MEDIUMTEXT,
    pw_ids MEDIUMTEXT,
    FOREIGN KEY (uniprot_id) REFERENCES uniprotkb(id),
    PRIMARY KEY (uniprot_id)
    );
        
/*A table that says if a table has been created or not. If yes it should be there.*/                              
CREATE TABLE goodTables (
    id int NOT NULL AUTO_INCREMENT,
    tableName varchar(50) not null,
    unique (tableName),
    PRIMARY KEY (ID) 
  );

/*A table containing info about the database*/                              
CREATE TABLE infoTable (
    dbName varchar (20),
    dbRelDate varchar (20),
    annocriptVer varchar (20),
    dbTablesVersion varchar (20)
  );
  
/*                        

create view ViewSeq as
      select uniprotId, seqId, OSType, description,  go_id
      FROM ( select u.id as uniprotId, u.seqId as seqId, g.go_id as go_id
                FROM  uniprotkb u LEFT JOIN mappingsTable g  ON u.id=g.uniprot_id
                ) 
            LEFT JOIN (select u.id as uniprotId2, u.seqId as seqId2, os.type as OSType, d.type as description
                FROM (uniprotkb u INNER JOIN OSType os ON u.OSId=os.id)
                        INNER JOIN  uniprotkbDescType d ON u.descId=d.id
            ) ON uniprotId=uniprotId2;
*/

/*ViewSeq is a view composed of the following two SELECT
---UNCOMMENT IT TO SEE THEM WORKING SEPARATELY---

*/

/*This is a composite view formed by two views (A, B) that permits to access the taxonomy, description and go ids of a 
given protein using its uniprot ID*/
create view ViewSeqA as
    select u.id as uniprotId, u.trId as trId, u.seqId as seqId, g.go_ids as go_ids
        FROM  uniprotkb u LEFT JOIN uniprotkbMappingsTable g  ON  g.uniprot_id = u.id;
           
create view ViewSeqB as
    select u.id as uniprotId2, os.type as OSType, d.type as description
	FROM (uniprotkb u INNER JOIN OSType os ON u.OSId=os.id)
                        INNER JOIN  uniprotkbDescType d ON u.descId=d.id;
                        
create view uniprotkbViewSeq as
    select uniprotId, trId, seqId, OSType, description, go_ids
        from ViewSeqA LEFT JOIN ViewSeqB ON uniprotId=uniprotId2;
        
        
/*This is a view that permits to access the enzyme ids of a 
given protein using its uniprot ID*/   
create view uniprotkbViewEnzyme as
    select u.seqId as seqId, m.enzyme_ids as enzyme_ids
        FROM  uniprotkbMappingsTable m LEFT JOIN uniprotkb u ON  m.uniprot_id = u.id;
  

/*This is a view that permits to access the pathways ids of a 
given protein using its uniprot ID*/   
create view uniprotkbViewPathways as
    select u.seqId as seqId, m.pw_ids as pw_ids
        FROM  uniprotkbMappingsTable m LEFT JOIN uniprotkb u ON  m.uniprot_id = u.id;
  
