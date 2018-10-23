#!/usr/bin/perl

#Annocript - A complete tool for transcriptomes annotation
#Copyright (C) <2018>  <Francesco Musacchia>

#This program is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.

#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.

#You should have received a copy of the GNU General Public License
#along with this program.  If not, see <http://www.gnu.org/licenses/>.
    
    
package LIB::mysql_db_management;
## mysql_db_management.pm
#Author: Francesco Musacchia  (2018)
#Permits the management of the database for the pipeline

#BEGIN { print (("  " x $main::x++) . "Beginning dbmanagement compile\n") }#DEBUGCODE

BEGIN
{
require Exporter;
use vars qw(@ISA @EXPORT);
@ISA = qw(Exporter);
@EXPORT_OK = qw( check_DB_user_permissions db_present destroy_db drop_database) ;
}
use strict;
use warnings;

use DBI;#DB library

use LIB::programs_management qw(log_and_exit print_and_log );

=head2 drop_database

 Title : drop_database
 Usage : drop_database( -file => 'dna2pepout_2out',
					);

 Function: creates a gff3 file from the portrait output and an hash with same information. 
          Example of Portrait output: >comp165_c0_seq1 :1: 0.948593 0.0514073

 Returns : a gff3 file

=cut
sub drop_database{
  my $db = shift;
  my $platform = shift;
  my $host = shift;
  my $port = shift;
  my $mySqlUser = shift;
  my $mySqlPass = shift;
  
  my $correctlyDestroyed = 0;
    
  #Connect to the information_schema db (to connect to some db)
  my $dsn = "dbi:$platform:information_schema:$host:$port";
  my $dbh = DBI->connect($dsn, $mySqlUser, $mySqlPass, { RaiseError => 1 } )
              or die "Unable to connect: $DBI::errstr\n";
  
  #Use the func subroutine to drop the database
  my $rc = $dbh->func("dropdb", $db, $host, $mySqlUser, $mySqlPass, 'admin');
  
  #Errors management   
  if ($dbh->err()) { 
    warn "$DBI::errstr. Unable to drop database $db\n"; 
  }
  else{ 
    $correctlyDestroyed = 1;
    print "$db correctly destroyed!\n";
  }
    
  $dbh->disconnect();
    
  return $correctlyDestroyed;
}


=head2 destroy_db

 Title   : destroy_db
 Usage   : destroy_db( -database => 'name of the database,
                               );

 Function:  destroy a database with all the tables

 Returns : nothing

=cut
sub destroy_db{
    my $database = shift;
    my $platform = shift;
    my $host = shift;
    my $port = shift;
    my $mySqlUser = shift;
    my $mySqlPass = shift;
    
    
    # DATA SOURCE NAME
    my $dsn = "dbi:$platform:information_schema:$host:$port";
    my $dbh = DBI->connect($dsn, $mySqlUser, $mySqlPass, { RaiseError => 1 } )
              or die "ERROR [$?]: unable to connect: $DBI::errstr\n";
    
    #First the database with that name will be erased
    $dbh->do("DROP DATABASE IF EXISTS $database") or die "ERROR: unable to cancel database $database! ?"; 
}



=head2 db_present

 Title   : db_present
 Usage   : db_present( -dbToSearch => 'name of the database);

 Function:  checks if a database exists.

 Returns : 1 if the database exists, 0 otherwise

=cut
sub db_present {
  my $dbToSearch = shift;
  my $platform = shift;
  my $host = shift;
  my $port = shift;
  my $mySqlUser = shift;
  my $mySqlPass = shift;
  
  my $dsn = "dbi:$platform:information_schema:$host:$port";
              
  my $dbh = DBI->connect($dsn, $mySqlUser, $mySqlPass, { RaiseError => 1} );
  my $databases = $dbh->selectcol_arrayref("show databases like '".$dbToSearch."'");
  
  my $numKeys = @$databases;
  

  return $numKeys;
}

=head2 check_DB_user_permissions

 Title   : check_DB_user_permissions
 Usage   : check_DB_user_permissions(  - local folder = the folder where actually is the PROGRAMS_OUT_PARSER script
                               );

 Function:  this function accesses the database with the user and password given and says if the user has sufficient privileges to 
        create a database
 
 Returns : 1 if check succeed -1 otherwise

=cut
sub check_DB_user_permissions{
	my $platform = shift;
	my $host = shift;
	my $port = shift;
	my $mySqlUser = shift;
	my $mySqlPass = shift;
  
  
  my $retVal = -1;
  my $privilegesType = 'GRANT ALL PRIVILEGES ON';
  
  print "Checking your privileges for MySQL...\n";
  # PERL DBI CONNECT
  my $dsn = "dbi:$platform:information_schema:$host:$port";
  my $dbh = DBI->connect($dsn, $mySqlUser, $mySqlPass, { RaiseError => 1 } ) 
    or die "Unable to connect (Check mySqlUser and mySqlPass in the config file): $DBI::errstr\n";;
  
  #Query to show grants for the user
  my $query = "SHOW GRANTS FOR '$mySqlUser'\@'$host'";
  my $sth = $dbh->prepare($query);
  $sth->execute();
  my  @data = $sth->fetchrow_array();
  
  my $privileges = $data[0];
  #print( "$privileges;\n\n");
   
  if($privileges =~ /$privilegesType/){
    $retVal = 1;
    print "DONE!\n";
  }
}

1;
