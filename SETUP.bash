#!/bin/bash

param="no"
check="OCCUPIED"
port=0000
my_array=()
RD='\033[0;31m'
GN='\033[1;32m'
BN='\033[0;33m'
BE='\033[0;34m'
PE='\033[0;35m'
NC='\033[0m' # No Color

printf "\n${BE}#=======================================================#${NC} \n\n"
printf "\n${BE}-----! Welcome to Docker Compose yaml Setup script !-----${NC} \n\n"
printf "\n${BE}#=======================================================#${NC} \n\n"


while [[ "${param,,}" == *"n"* ]]
do
	printf "Please provide the project name(container name) :\n(Use '-' as separator)\n\n"
	read -p 'Project Name : ' pronmvar
	# -- Pre process -- #
	pronm=${pronmvar//[^[:alnum:]-]/}
	pro=${pronm,,}
	if [[ -z $pro ]]
	then
		printf "\n${RD}Please provide a valid input${NC}\n"
                exit
	fi
	
	#printf "\nWhich directory would you like to the map to the project Directory of the container : [W]orking directory or [C]ustom or [N]one ?(w/c/n) :"
	printf "\nMount the working directory to the project Directory of the container ?(y/n) :"
	read pdir

	if [[ "${pdir,,}" == *"y"* ]]
	then
		spro="."
		my_array+=( "$spro" )
		my_array+=( "/home/rstudio/project" )
	elif [[ "${pdir,,}" == *"n"* ]]
	then
		spro=""
	else
		printf "\n${RD}Please provide a valid input${NC}\n"
		exit
	fi

	echo "Mount the Default R .cache {\$HOME/.cache} or Custom or None ?(y/c/n) :"
	read dir

	if [[ "${dir,,}" == *"y"* ]]
	then
		scache="\$HOME/.cache"
		my_array+=( "$scache" )
		my_array+=( "/home/rstudio/.cache" )
	elif [[ "${dir,,}" == *"c"* ]]
	then
		read -p '.cache directory : ' scache
		my_array+=( "$scache" )
		my_array+=( "/home/rstudio/.cache" )
	elif [[ "${dir,,}" == *"n"* ]]
	then
		scache=""
	else
		printf "\n${RD}Please provide a valid input${NC}\n"
		exit
	fi

	echo "How many additional directories would you like to mount? : "
	read dirnum

	while [[ $dirnum > 0 ]]; do
	echo "Please provide the Source{Host System} and Destination{docker image} directories with Full Path"
		read -p 'Source : ' scdir
		read -p 'Destination : ' dsdir
		my_array+=( "$scdir" )
		my_array+=( "$dsdir" )
		((dirnum--))
	done


	while [[ "$check" == "OCCUPIED" ]]
	do
		lowerRange=50000   # inclusive
		upperRange=60000   # exclusive
		port=$(( RANDOM * ( upperRange - lowerRange) / 32767 + lowerRange ))
		check=$(nc -z 127.0.0.1 $port && echo "OCCUPIED" || echo "FREE")
	done

	printf "\nWould you like keep the Default password{1rstudio} or Change or None ?(y/c/n) :"
	read pasr

	if [[ "${pasr,,}" == *"y"* ]]
	then
		spas="1rstudio"
	elif [[ "${pasr,,}" == *"c"* ]]
	then
		read -p 'Please provide a password : ' spas
	elif [[ "${pasr,,}" == *"n"* ]]
	then
		spas=""
	else
		printf "\n${RD}Please provide a valid input${NC}\n"
		exit
	fi
	# ----------------- #
	printf "\n${PE}Please Review the Parameters below${NC}\n\n"
	printf "\n${BN}Project Name :${NC} $pro"
	printf "\n${BN}Mounted R .cache :${NC} $scache"
	arraylength=${#my_array[@]}
	for (( i=0; i<${arraylength}; i+=2 ));
	do
		printf "\n${BN}Mounted Directories :${NC} ${GN}Source${NC} = ${my_array[$i]} || ${BE}Destination${NC} = ${my_array[$i+1]}"
	done

	printf "\n${BN}Mounted Port :${NC} $port $check"
	printf "\n${BN}USERNAME :${NC} $USER"
	printf "\n${BN}USERID :${NC} $(id -u)"
	printf "\n${BN}GROUPID :${NC} $(id -g)"

	printf "\n\n${GN}Do you accept the parameters?${NC} (Yes/No/Discard)"
	read param
done


if [[ "${param,,}" == *"y"* ]]
then
	printf "services:
  rstudio:
    image: rstudiorenvtidy:4.4.3
    build:
     context: ./docker
    container_name: $pro
    ports:
      - \"$port:8787\"
    volumes:
      - type: \"bind\"
        source: \"./config\"
        target: \"/home/rstudio/.config/rstudio\"" > compose.yml
	for (( i=0; i<${arraylength}; i+=2 ));
        do
                printf "
      - type: \"bind\"
        source: \"${my_array[$i]}\"
        target: \"${my_array[$i+1]}\"" >> compose.yml
        done
	printf "
    environment:
      - USERNAME=$USER
      - USERID=$(id -u)
      - GROUPID=$(id -g)
      - USER=rstudio
      - PASSWORD=$spas" >> compose.yml
      printf "\n${RD}Would you like me to start the container?${NC}(y/n)"
      read run

        if [[ "${run,,}" == *"y"* ]]
	then
		docker compose -f ./compose.yml up -d
		printf "\nPlease open browser and browse to ${GN}http://localhost:$port/${NC}\n\nRstudio USER     : rstudio\nRstudio PASSWORD : $spas\n\n"
	fi


fi

#docker compose -f ./compose.yml up -d
