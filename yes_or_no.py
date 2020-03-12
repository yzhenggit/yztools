import sys 

def yes_or_no(question):
    while True:
        sys.stdout.write(question + '[y/n]\n')
        answer = input()
        if answer == 'y':
            return True
        elif answer == 'n':
            return False
        else:
            sys.stdout.write("Please respond with y or n \n")

