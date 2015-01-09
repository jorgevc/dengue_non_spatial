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
