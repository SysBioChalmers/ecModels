pipeline {
  agent any
  stages {
    stage('Run') {
//      when {
//        branch 'develop'
//      }
      steps {
        sh '''
          python3 run.py
          '''
      }
    }
  }
}
