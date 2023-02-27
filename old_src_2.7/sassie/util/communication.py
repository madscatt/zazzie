'''
communication.py

contains methods to communicate with GUIs

'''

import json, sys, StringIO, socket

def tcpquestion(json_variables, question, timeout=300, buffersize=65536):
    '''
    method to communicate to SASSIE-web via GenApp QRM TCP protocol

    # json_variables is the input json object
    # question is either a dict or json string
    # timeout is in seconds
    # buffersize is for the answer, so if you expect larger than 64k of total json string size, use a larger number


    '''

    if '_tcphost' not in json_variables or '_tcpport' not in json_variables:
        return json.loads('{"error":"no tcp support"}')

    # build question

    msg = {}
    msg['_uuid'] = json_variables['_uuid']
    msg['timeout'] = timeout
    if isinstance( question, basestring ):
        try:
            msg['_question'] = json.loads(question)
        except ValueError:
            return json.loads('{"error":"json question decoding error"}')
    elif isinstance( question, dict ):
        msg['_question'] = question
    else:
        return json.loads('{"error":"question must be a json string or dict"}')

    msgj = json.dumps(msg)
    # a newline is also required when sending a question
    msgj += '\n'

    # send question
    s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    s.connect((json_variables['_tcphost'],json_variables['_tcpport']))
    s.send(msgj)

    # receive answer

    data = s.recv(buffersize)
    s.close()
    return json.loads(data)

if __name__=='__main__':

        '''
        main method to test tcpquestion

        '''

        argv_io_string = StringIO.StringIO(sys.argv[1])
        json_variables = json.load(argv_io_string)

        mass = float(json_variables['m'])
        speed_of_light = float(json_variables['c'])

        import mass_energy

        output = {}

        # ask silly question if possible

        myquestion = '''
{
    "id" : "q1"
    ,"title" : "are you sure?"
    ,"text" : "<p>header text.</p><hr>"
    ,"fields" : [
        {
            "id" : "l1"
            ,"type" : "label"
            ,"label" : "<center>this is label text</center>"
        }
        ,{
            "id" : "t1"
            ,"type" : "text"
            ,"label" : "tell me your name:"
        }
        ,{
            "id" : "cb1"
            ,"type" : "checkbox"
            ,"label" : "are you sure about the speed of light?"
        }
    ]
}
'''.strip()
        answer = tcpquestion( json_variables, myquestion );

# alternatively, one could pass a dict
#        mydict = json.loads(myquestion)
#        ... manipulations of mydict ...
#        answer = tcpquestion( json_variables, mydict );

        if 'error' in answer:
                # handle error, could be timeout, no tcp support, malformed question json or other
                output['question_error'] = answer['error']

        if answer['_response']['button'] == 'cancel':
                # the user pressed cancel, handle it
                output['user_input'] = 'user has canceled'

        output['users_name'] = answer['_response']['t1']

        output['e'] = mass_energy.einstein(mass,speed_of_light)
        output['answer'] = json.dumps(answer)

        print json.dumps(output)





