/*
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 * MA 02110-1301, USA.
 * 
 */

#include "stdio.h"
#include "conn_mysql.h"


void finish_with_error(MYSQL *con)
{
  fprintf(stderr, "%s\n", mysql_error(con));
  mysql_close(con);
  return;        
}

MYSQL* connect_db(char *host,char *user, char *passw, char *db)
{
MYSQL *con = mysql_init(NULL);

if (con == NULL) 
  {
      fprintf(stderr, "%s\n", mysql_error(con));
  }  

  if (mysql_real_connect(con, host, user, passw, 
          db, 0, NULL, 0) == NULL) 
  {
      finish_with_error(con);
  }
 
return con;
}
  
int insert_into_db(MYSQL *con, char *table,char *values)
{ 
	
	char query[200];
	sprintf(query,"INSERT INTO %s VALUES(%s)",table,values);
  if (mysql_query(con, query)) {
      finish_with_error(con);
  }
  
return mysql_insert_id(con);
}
