<?php


/*
|--------------------------------------------------------------------------
| Configure your contact form
|--------------------------------------------------------------------------
|
| Set value of '$reciever' to email address that want to receive inquiries.
| Also, '$default_subject' is the subject that you'll see in your inbox.
|
*/
$reciever = "info@yourdomain.com";
$default_subject = "[Yourdomain.com] - Inquery from " . strip_tags($_POST['name']);



/*
|--------------------------------------------------------------------------
| Sending email
|--------------------------------------------------------------------------
|
| This part of code is responsible to send the email. So you don't need to
| change anything here. Just make sure you have a mailserver on your host.
|
*/

if ( isset($_POST['email']) && isset($_POST['message']) && filter_var($_POST['email'], FILTER_VALIDATE_EMAIL) ) {

  // detect & prevent header injections
  $test = "/(content-type|bcc:|cc:|to:)/i";
  foreach ( $_POST as $key => $val ) {
    if ( preg_match( $test, $val ) ) {
      exit;
    }
  }

  $subject = $_POST['subject'];
  if ($subject == "") {
    $subject = $default_subject;
  }

  $msg = nl2br($_POST['message']);

  $headers = "From: " . strip_tags($_POST['name']) . '<'. strip_tags($_POST['email']) .'>' . "\r\n";
  $headers .= "Reply-To: ". strip_tags($_POST['email']) . "\r\n";
  $headers .= "MIME-Version: 1.0\r\n";
  $headers .= "Content-Type: text/html; charset=ISO-8859-1\r\n";
  
  mail( $reciever, $subject , $msg, $headers );
}
?>