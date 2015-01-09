#include <mysql.h>

void finish_with_error(MYSQL *con);

MYSQL* connect_db(char *host,char *user, char *passw, char *db);

int insert_into_db(MYSQL *con, char *table,char *values);
