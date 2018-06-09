
while getopts p:f:o: option
do
	case "${option}"
		in
		p) PRIVKEY=${OPTARG};;
		f) FILE=$OPTARG;;
		o) OUT=$OPTARG;;
	esac
done


#Entlüsseln des Key-Files "key.bin" mit dem Faktorisieren Private Key
echo "Decrypting the Key-File with the factorized RSA-Private-Key"
openssl rsautl -decrypt -inkey $PRIVKEY -in key.bin.enc -out key.bin

# Entschlüsseln der AES verschlüsselten Datei mithilfe des entschlüsselten Key-Files:
echo "Decrypting the File '' with the decrypted Key-File"
openssl enc -d -aes-256-cbc -in $FILE -out $OUT -pass file:./key.bin
