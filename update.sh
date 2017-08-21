system('scp -r * ps562@ws190.gurdon.private.cam.ac.uk:/var/shiny-server/www/DEview')

system('rsync -avp server.R ui.R ps562@ws190.gurdon.private.cam.ac.uk:/var/shiny-server/www/DEview')


system('ssh ps562@ws190 ls /var/shiny-server/www/DEview')